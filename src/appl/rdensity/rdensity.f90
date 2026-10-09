! Shared entry point for rdensity and rdensity_mpi.
program rdensity
    use default_C, only: NDEF
    use debug_C, only: CUTOFF
    use density_parallel
    use density_calculation
    use getyn_I
    use setdbg_I
    use setmc_I
    use setcon_I
    use setcsla_I
    use getsmd_I
    use getmixblock_I
    use factt_I
    implicit none
    character(24) :: name
    character(256) :: arg
    integer :: nci, ordering, status, i, ios, length, ncore
    logical :: yes, exists, interactive

    call parallel_start
    name = ''
    nci = 1
    ndef = 0
    ordering = 2
    status = 0
    CUTOFF = 1.0D-10
    interactive = command_argument_count() == 0
    if (rank == 0) then
        if (interactive) then
            write(*,*) 'Default settings?'
            yes = getyn()
            if (.not.yes) ndef = 1
            if (workers > 1 .and. ndef /= 0) then
                write(*,*) 'MPI runs support default settings only, as in RHFS/RIS.'
                status = 2
            else
                write(*,*) 'Name of state'
                read(*,'(a)',iostat=ios) arg
                if (ios /= 0) status = 2
                call accept_name(arg)
                write(*,*) 'Mixing coefficients from a CI calc.?'
                yes = getyn()
                if (yes) nci = 0
            endif
        else
            do i = 1,command_argument_count()
                call get_command_argument(i, arg, length=length, status=ios)
                if (ios /= 0 .or. length > len(arg)) then
                    status = 2
                    exit
                endif
                select case(trim(arg))
                case('--ci')
                    nci = 0
                case('--nonci')
                    nci = 1
                case('--order=1')
                    ordering = 1
                case('--order=2')
                    ordering = 2
                case('--help','-h')
                    status = 1
                    exit
                case default
                    if (i == 1 .and. arg(1:1) /= '-') then
                        call accept_name(arg)
                    else
                        write(*,*) 'Unknown argument: ', trim(arg)
                        status = 2
                    endif
                end select
            enddo
        endif
        if (status == 0 .and. len_trim(name) == 0) status = 2
        if (status /= 0) write(*,*) 'Usage: rdensity[_mpi] name [--ci|--nonci] [--order=1|--order=2]'
    endif
    call broadcast_options(name, nci, ndef, ordering, status)
    if (status == 1) then
        call parallel_finish
        stop
    endif
    if (status /= 0) call density_fail('Invalid command line or interactive input')

    call require_file('isodata')
    call require_file(trim(name)//'.c')
    call require_file(trim(name)//'.w')
    if (nci == 0) then
        call require_file(trim(name)//'.cm')
        arg = trim(name)//'.cd'
    else
        call require_file(trim(name)//'.m')
        arg = trim(name)//'.d'
    endif
    if (rank == 0) then
        inquire(file=trim(arg), exist=exists)
        if (exists) call density_fail('Density output already exists: '//trim(arg))
    endif
    call setdbg
    call setmc
    call setcon
    ! Read-only inputs are shared; angular and orbital state is rank-local.
    call setcsla(name, ncore)
    call getsmd(name)
    call getmixblock(name, nci)
    call factt
    if (rank == 0 .and. interactive) then
        write(*,*) 'How do you want to order your eigenvectors?'
        write(*,*) '1) Dominant component; 2) Decreasing occupation (preserving pure orbitals)'
        read(*,*,iostat=ios) ordering
        if (ios /= 0) status = 2
        if (ordering /= 1 .and. ordering /= 2) status = 2
    endif
    call broadcast_options(name, nci, ndef, ordering, status)
    if (status /= 0) call density_fail('Orbital ordering must be 1 or 2')
    call calculate_density(name, nci, ordering)
    call parallel_finish
contains
    subroutine accept_name(text)
        character(*), intent(in) :: text
        ! Legacy readers find the trailing blank in CHARACTER(24).
        if (len_trim(text) < 1 .or. len_trim(text) > 23) then
            status = 2
        else if (index(trim(text),' ') /= 0) then
            status = 2
        else
            name = trim(text)
        endif
    end subroutine

    subroutine require_file(filename)
        character(*), intent(in) :: filename
        logical :: found
        inquire(file=filename, exist=found)
        if (.not.found) call density_fail('Missing input file: '//filename)
    end subroutine
end program rdensity
