! Based on GRASP src/tool/rlevels.f90 (Xinghong He, 1998; G. Gaigalas, 2011).
program rlevels
    ! Module dependencies
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: error_unit, int64, iostat_end, output_unit, real64
    use calculation_state_files, only: read_state_entry, resolve_mixing_file
    implicit none (type, external)

    real(real64), parameter :: rydberg_wavenumber_cm = 109737.31568508_real64
    real(real64), parameter :: hartree_to_wavenumber_cm = 2.0_real64*rydberg_wavenumber_cm
    real(real64), parameter :: label_energy_tolerance = 1.0e-8_real64
    integer, parameter :: label_column = 37, label_width = 164 ! jj2lsj_omp .lsj.lbl layout

    ! Derived types
    type :: level
        ! Components
        integer :: root_position_in_symmetry_block
        real(real64) :: energy
        character(len=12) :: j_text
        character :: parity
        character(len=label_width) :: configuration = ''
    end type level

    ! Derived types
    type :: symmetry_block
        ! Components
        integer :: first_level_index = 0, block_level_count = 0
    end type symmetry_block

    type(level), allocatable :: levels(:)
    character(len=:), allocatable :: path
    integer :: argument_index, path_length, io_status, level_count
    logical :: show_configuration

    level_count = 0
    show_configuration = .false.
    allocate(levels(0))
    select case (command_argument_count())
    case (0)
        write (error_unit, '(A)') 'Usage: rlevels file1.mix file2.cimix ...', &
            'Enter one state name or mixing file per line (blank line or EOF to finish).'
        do
            call read_state_entry(path, io_status)
            if (io_status == iostat_end .and. len(path) == 0) exit
            if (io_status /= 0 .and. io_status /= iostat_end) call fail('standard input', 'cannot read state name or mixing file')
            path = trim(adjustl(path))
            if (len(path) == 0) exit
            call read_mixing_file(resolve_mixing_file(path, interactive=.true.), levels, level_count, show_configuration)
            if (io_status == iostat_end) exit
        end do
    case default
        do argument_index = 1, command_argument_count()
            call get_command_argument(argument_index, length=path_length, status=io_status)
            if (io_status /= 0 .or. path_length == 0) call fail('command line', 'expected a mixing file name')
            allocate(character(len=path_length) :: path)
            call get_command_argument(argument_index, value=path, status=io_status)
            if (io_status /= 0) call fail('command line', 'cannot read mixing file name')
            call read_mixing_file(resolve_mixing_file(path, interactive=.false.), levels, level_count, show_configuration)
            deallocate(path)
        end do
    end select
    if (level_count == 0) call fail('input', 'no energy levels found')
    call print_levels(levels(:level_count), show_configuration)

contains

    subroutine read_mixing_file(name, entries, used, with_labels)
        ! Arguments
        character(len=*), intent(in) :: name
        type(level), allocatable, intent(inout) :: entries(:)
        integer, intent(inout) :: used
        logical, intent(inout) :: with_labels

        ! Local declarations
        type(level), allocatable :: grown(:)
        type(symmetry_block), allocatable :: blocks(:)
        character(len=6) :: magic
        character(len=12) :: j_text
        character :: parity
        integer :: unit, status, electron_count, csf_count, orbital_count, file_level_count, mixing_coefficient_count, &
            & symmetry_block_count
        integer :: symmetry_block_index, last_block, block_csf_count, block_level_count, two_j_plus_one, parity_code, first, &
            & last, file_end, suffix
        integer(int64) :: coefficient_index
        real(real64) :: hamiltonian_energy_shift, discarded_coefficient
        logical :: mcdhf

        suffix = index(name, '.', back=.true.)
        mcdhf = .false.
        if (suffix > 0) mcdhf = name(suffix:) == '.mix' .or. name(suffix:) == '.m'
        open (newunit=unit, file=name, form='unformatted', status='old', action='read', iostat=status)
        if (status /= 0) call fail(name, 'cannot open mixing file')
        read (unit, iostat=status) magic
        if (status /= 0) call fail(name, 'cannot read G92MIX signature')
        if (magic /= 'G92MIX') call fail(name, 'expected G92MIX signature')
        read (unit, iostat=status) electron_count, csf_count, orbital_count, file_level_count, mixing_coefficient_count, &
            & symmetry_block_count
        if (status /= 0) call fail(name, 'cannot read mixing header')
        if (file_level_count < 0 .or. symmetry_block_count < 0 .or. file_level_count > huge(used) - used) &
            call fail(name, 'invalid mixing header counts')

        file_end = used + file_level_count
        if (file_end > size(entries)) then
            allocate(grown(max(file_end, size(entries) + min(size(entries), huge(used) - size(entries)))), stat=status)
            if (status /= 0) call fail(name, 'cannot allocate energy levels')
            grown(:used) = entries(:used)
            call move_alloc(grown, entries)
        end if
        allocate(blocks(symmetry_block_count), stat=status)
        if (status /= 0) call fail(name, 'cannot allocate symmetry blocks')
        last = used
        last_block = 0
        do
            read (unit, iostat=status) symmetry_block_index, block_csf_count, block_level_count, two_j_plus_one, parity_code
            if (status == iostat_end) exit
            if (status /= 0) call fail(name, 'cannot read symmetry block header')
            ! RMCDHF omits blocks with no requested roots. RCI writes their headers.
            if (symmetry_block_index <= last_block .or. symmetry_block_index > symmetry_block_count .or. (.not. mcdhf .and. &
                & symmetry_block_index /= last_block + 1)) &
                call fail(name, 'unexpected symmetry block number')
            last_block = symmetry_block_index
            if (block_level_count < 0 .or. block_level_count > block_csf_count .or. block_level_count > file_end - last) &
                call fail(name, 'invalid root count in symmetry block')
            if (block_level_count == 0) cycle
            if (two_j_plus_one < 1 .or. (parity_code /= 1 .and. parity_code /= -1)) &
                call fail(name, 'invalid J or parity in symmetry block')
            j_text = format_j(two_j_plus_one - 1)
            parity = merge('+', '-', parity_code == 1)
            first = last + 1
            last = last + block_level_count
            blocks(symmetry_block_index) = symmetry_block(first, block_level_count)
            read (unit, iostat=status) entries(first:last)%root_position_in_symmetry_block
            if (status /= 0) call fail(name, 'cannot read root positions')
            if (any(entries(first:last)%root_position_in_symmetry_block < 1)) call fail(name, 'invalid root position')
            read (unit, iostat=status) hamiltonian_energy_shift, entries(first:last)%energy
            if (status /= 0) call fail(name, 'cannot read block energies')
            if (.not. ieee_is_finite(hamiltonian_energy_shift) .or. .not. all(ieee_is_finite(entries(first:last)%energy))) &
                call fail(name, 'nonfinite block energy')
            ! Consume the coefficient record without storing the eigenvectors.
            read (unit, iostat=status) &
                (discarded_coefficient, coefficient_index=1_int64, int(block_csf_count, int64)*block_level_count)
            if (status /= 0) call fail(name, 'cannot read mixing coefficient record')
            entries(first:last)%energy = entries(first:last)%energy + hamiltonian_energy_shift
            entries(first:last)%j_text = j_text
            entries(first:last)%parity = parity
        end do
        close (unit)
        if (last /= file_end .or. (.not. mcdhf .and. last_block /= symmetry_block_count)) &
            call fail(name, 'mixing records disagree with header counts')

        if (suffix > 0 .and. file_level_count > 0) then
            select case (name(suffix:))
            case ('.mix', '.cimix', '.m', '.cm')
                call read_labels(name(:suffix-1)//'.lsj.lbl', entries, blocks, with_labels)
            end select
        end if
        used = file_end
        write (output_unit, *) 'nblock = ', symmetry_block_count, '  ncftot = ', csf_count, '  nw = ', orbital_count, &
            & '  nelec = ', electron_count
    end subroutine read_mixing_file

    function format_j(two_j) result(j_text)
        ! Arguments
        integer, intent(in) :: two_j
        character(len=12) :: j_text

        if (mod(two_j, 2) == 0) then
            write (j_text, '(I0)') two_j/2
        else
            write (j_text, '(I0,A)') two_j, '/2'
        end if
    end function format_j

    subroutine read_labels(name, entries, blocks, with_labels)
        ! Arguments
        character(len=*), intent(in) :: name
        type(level), intent(inout) :: entries(:)
        type(symmetry_block), intent(in) :: blocks(:)
        logical, intent(inout) :: with_labels

        ! Local declarations
        character(len=label_column+label_width-1) :: line
        character(len=4) :: j_text
        character :: parity
        real(real64) :: energy
        integer :: unit, status, ordinal, block, index_in_table
        logical :: exists

        inquire (file=name, exist=exists, iostat=status)
        if (status /= 0) call fail(name, 'cannot check for LSJ labels')
        if (.not. exists) return
        open (newunit=unit, file=name, form='formatted', status='old', action='read', iostat=status)
        if (status /= 0) call fail(name, 'cannot open LSJ labels')
        read (unit, '(A)', iostat=status) line
        if (status /= 0) call fail(name, 'cannot read LSJ label header')
        if (index(line, 'Energy Total') == 0) call fail(name, 'invalid LSJ label header')
        with_labels = .true.
        do
            read (unit, '(A)', iostat=status) line
            if (status == iostat_end) exit
            if (status /= 0) call fail(name, 'cannot read LSJ label record')
            ! An I3 position can overflow in the producer. That row has no usable ordinal.
            if (line(:3) == '' .or. line(:3) == '***') cycle
            read (line, '(I3,1X,A4,5X,A1,8X,F16.9)', iostat=status) ordinal, j_text, parity, energy
            if (status /= 0) call fail(name, 'invalid LSJ level fields')
            if (ordinal < 1 .or. .not. ieee_is_finite(energy)) call fail(name, 'invalid LSJ position or energy')
            read (unit, '(A)', iostat=status) line
            if (status /= 0) call fail(name, 'cannot read leading LSJ component')
            if (line(:7) /= '' .or. len_trim(line) < label_column) call fail(name, 'missing leading LSJ configuration')

            ! jj2lsj prints the ordinal in the block, which can differ from the root position in G92MIX.
            do block = 1, size(blocks)
                if (ordinal > blocks(block)%block_level_count) cycle
                index_in_table = blocks(block)%first_level_index + ordinal - 1
                associate (entry => entries(index_in_table))
                    if (entry%j_text /= adjustl(j_text) .or. entry%parity /= parity) cycle
                    if (abs(entry%energy - energy) >= label_energy_tolerance) cycle
                    entry%configuration = line(label_column:)
                end associate
                exit
            end do
        end do
        close (unit)
    end subroutine read_labels

    recursive subroutine sort_by_energy(entries, order, work)
        ! Arguments
        type(level), intent(in) :: entries(:)
        integer, intent(inout) :: order(:), work(:)

        ! Local declarations
        integer :: half, left, right, k

        half = size(order)/2
        if (half == 0) return
        call sort_by_energy(entries, order(:half), work(:half))
        call sort_by_energy(entries, order(half+1:), work(half+1:))
        left = 1
        right = half + 1
        do k = 1, size(order)
            if (left > half .or. right > size(order)) exit
            if (entries(order(left))%energy <= entries(order(right))%energy) then
                work(k) = order(left)
                left = left + 1
            else
                work(k) = order(right)
                right = right + 1
            end if
        end do
        if (left <= half) then
            work(k:) = order(left:half)
        else
            work(k:) = order(right:)
        end if
        order = work
    end subroutine sort_by_energy

    subroutine print_levels(entries, with_labels)
        ! Arguments
        type(level), intent(in) :: entries(:)
        logical, intent(in) :: with_labels

        ! Local declarations
        character(len=:), allocatable :: rule, heading
        character(len=4) :: legacy_j
        integer, allocatable :: order(:), work(:)
        integer :: row, legacy_j_width
        real(real64) :: ground_energy_hartree, previous_level_energy_hartree, excitation, splitting

        allocate(order(size(entries)), work(size(entries)))
        order = [(row, row=1,size(entries))]
        call sort_by_energy(entries, order, work)

        write (output_unit, '(A)') '', ' Energy levels for ...'
        write (output_unit, '(A,F14.5)') ' Rydberg constant is ', rydberg_wavenumber_cm
        heading = ' No Pos  J Parity Energy Total    Levels     Splitting'
        if (with_labels) then
            heading = heading//'     Configuration'
            rule = repeat('-', 90)
            write (output_unit, '(A)') ' Splitting is the energy difference with the lower neighbor'
        else
            rule = repeat('-', 73)
            write (output_unit, '(A)') ' No - Serial number of the state; Pos - Position of the state within the ', &
                ' J/P block; Splitting is the energy difference with the lower neighbor'
        end if
        write (output_unit, '(A)') rule, heading, '                      (a.u.)      (cm^-1)     (cm^-1)', rule
        ground_energy_hartree = entries(order(1))%energy
        previous_level_energy_hartree = ground_energy_hartree
        do row = 1, size(order)
            associate (entry => entries(order(row)))
                excitation = (entry%energy - ground_energy_hartree)*hartree_to_wavenumber_cm
                splitting = (entry%energy - previous_level_energy_hartree)*hartree_to_wavenumber_cm
                legacy_j_width = merge(4, 3, index(entry%j_text, '/') > 0)
                if (row <= 999 .and. entry%root_position_in_symmetry_block <= 99 .and. &
                    & len_trim(entry%j_text) <= legacy_j_width .and. abs(entry%energy) <= 99999.0_real64 .and. &
                    & excitation < 1.0e7_real64) then
                    legacy_j = adjustr(entry%j_text(:4))
                    if (scan(legacy_j, '/') == 0) legacy_j = legacy_j(2:)//' '
                    write (output_unit, '(2I3,1X,A4,1X,A1,2X,F14.7,2F12.2,2X,A)') &
                        row, entry%root_position_in_symmetry_block, legacy_j, entry%parity, entry%energy, excitation, splitting, &
                            & trim(entry%configuration)
                else
                    write (output_unit, '(I0,1X,I0,1X,A,1X,A1,1X,F0.7,1X,F0.2,1X,F0.2,2X,A)') &
                        row, entry%root_position_in_symmetry_block, trim(entry%j_text), entry%parity, entry%energy, &
                        excitation, splitting, trim(entry%configuration)
                end if
                previous_level_energy_hartree = entry%energy
            end associate
        end do
        write (output_unit, '(A)') rule
    end subroutine print_levels

    subroutine fail(name, message)
        ! Arguments
        character(len=*), intent(in) :: name, message

        write (error_unit, '(A)') 'rlevels: '//name//': '//message
        error stop 1, quiet=.true.
    end subroutine fail
end program rlevels
