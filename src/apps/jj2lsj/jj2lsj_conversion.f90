! Transform jj-coupled CSF expansions into LSJ expansions and write spectroscopic labels.
!     This module contains the procedures which are used to performe
!     the jj-LSJ transformation.
!
! Written by G. Gaigalas,
! Vilnius                                  last update: Jan 2017
! Modified by Chongyang Chen. 2023
!     Some minor fixes labelled by CYC-2023
! Modified by G. Gaigalas. 2024
! Modified by G. Gaigalas. 2025
module jj2lsj_conversion
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use jj2lsj_profile, only: profile_emit, profile_file_size, profile_finish, profile_start, profile_time
    use jj2lsj_postprocessing, only: ls_ranking, post_worker_count, rank_ls_components
    implicit none (type, external)

    ! Module interface
    private :: real64
    ! R o u t i n e s
    public  :: transform_asf_to_ls
    ! Expands an atomic state function, which is represented
    ! in a jj-coupling CSF basis into a basis of LS-coupling CSF.
    public  :: ls_jj_recoupling_coefficient
    ! Returns the value of the LS-jj transformation matrix
    ! for a given set of quantum numbers.
    private :: single_partner_ls_jj_coefficient
    ! Returns the value of the LS-jj transformation matrix
    ! (l^2 LSJ| j_1 j_2 J).
    private :: lookup_ls_jj_subshell_coefficient
    ! Returns the value of the LS-jj transformation matrix
    ! form an array LS_jj_*_*.
    public  :: release_ls_expansion
    ! Dellocates the storage of asf_set_LS.
    public  :: format_ls_csf
    ! A spectroscopic notation of shell in LS coupling is return.
    public  :: format_ls_csf_label
    ! A spectroscopic notation of shell in LS coupling is return.
    public  :: select_ls_transition_labels
    public  :: enumerate_ls_subshell_terms
    ! This procedure return all allowed subshell terms
    ! (l, w, Q, L, S) for given l^N which must be 0, 1, 2 or 3.
    public  :: read_ls_conversion_options
    ! The input from the screen.
    public  :: parse_level_selection
    ! Attempts to interpret the serial level numbers from a
    ! string.
    public  :: run_ls_coupling_conversion
    ! Controls the transformation of atomic states from a jj-
    ! into a LS-coupled CSF basis.
    public  :: pack_ls_csf_label
    ! Encording all CSF in LS - coupling.
    public  :: write_jj_csf
    ! Prints all information about the single CSF in jj-coupling
    public  :: write_ls_csf
    ! Prints all information about the single CSF in LS-coupling
    public  :: write_ls_csf_list
    ! Prints all information about the single CSF scheme csf_set
    ! in LS-coupling
    public  :: build_ls_basis  ! Basis construction for the conversion interface.
    ! This subroutine fills up the variable asf_set_LS%csf_set_LS
    ! with data generated using the one from asf_set%csf_set
    !
    ! This subroutine contains the following internal routines:
    !  * subroutine setLS_action
    !    The subroutine defines the "action" of subroutine
    !    setLS_job_count: whether it counts
    !    the number of csfs_LS (asf_set_LS%csf_set_LS%novcsf)
    !    or fills the arrays of wave functions in LS coupling
    !    with asf_set_LS%csf_set_LS%csf(...) with
    !    the corresponding quantum nubers.
    !  * subroutine setLS_add_quantum_numbers
    !    The subroutine adds quantum numbers stored
    !    in temprorary arrays Li, Si, L_i, S_i, w, Q to
    !    the corresponding arrays of asf_set_LS%csf_set_LS%csf().
    !    private :: setLS_job_count
    !  * recursive subroutine setLS_job_count
    !    Recursive subroutine for the calculation of the
    !    number of csfs_LS and corresponding quantum numbers.
    !  * function setLS_equivalent_csfs
    !    This subroutine defines the "equivalency" of two csfs_jj
    !    in the sence of generation of csfs_LS
    !    number of csfs_LS and corresponding quantum numbers.
    public  :: ls_jj_csf_overlap
    ! Return the value of the transformation matrix
    ! from jj- to LS-coupling scheme in the case of any
    ! number of open shells.
    public  :: ls_jj_shell_recoupling_factor
    ! Return the value of main part of the transformation
    ! matrix from jj- to LS-coupling scheme in the
    ! case of any number of open shells.
    private :: assign_unique_lsj_labels

    ! Derived types
    !               Subroutine defines a unique labels for energy levels
    !
    !   D e f i n i t i o n   o f   A S F
    !   i n   L S -   C o u p l i n g
    type, public :: nonrelativistic_subshell_label
        ! Components
        integer :: n, l
    end type nonrelativistic_subshell_label

    ! Derived types
    type, public :: cs_function_ls
        ! Components
        integer :: two_j
        character(len=1)  :: parity
        integer, pointer :: occupation(:) => null()
        integer, pointer :: seniority(:) => null()
        integer, pointer :: w(:) => null()
        integer, pointer :: shell_two_l(:) => null()
        integer, pointer :: shell_two_s(:) => null()
        integer, pointer :: intermediate_two_l(:) => null()
        integer, pointer :: intermediate_two_s(:) => null()
    end type cs_function_ls

    ! Derived types
    type, public :: csf_basis_ls
        ! Components
        integer :: ls_csf_count         ! Number of CSF in the basis.
        integer :: shell_count      ! Number of (nonrelativistic) shells.
        integer :: closed_shell_count        ! Number of (closed) core shells.
        integer :: number_of_electrons
        type(nonrelativistic_subshell_label), pointer :: shell(:) => null()
        type(cs_function_ls), pointer :: csf(:) => null()
        type(parent_from_jj), pointer :: parent(:) => null()
    end type csf_basis_ls

    ! Derived types
    type, public :: as_function_ls
        ! Components
        integer          :: level_no
        integer          :: max_csf_index
        integer          :: totall, totals, two_j
        character(len=1) :: parity
        real(real64)     :: energy
        real(real64), pointer :: eigenvector(:) => null()
    end type as_function_ls

    ! Derived types
    type, public :: asf_basis_ls
        ! Components
        integer      :: asf_count           ! Number of considered ASF.
        real(real64) :: average_energy  ! Averaged energy of this set of ASF.
        type(as_function_ls), pointer :: asf(:) => null()
        type(csf_basis_ls) :: csf_set_ls
    end type asf_basis_ls

    ! Derived types
    !
    !   G l o b a l   V a r i a b l e s
    type, public :: parent_from_jj
        ! Components
        integer :: parent_minus
        integer :: parent_plus
    end type parent_from_jj

    ! Derived types
    type, public :: lsj_list
        ! Components
        integer::list_size                    !number items in a list
        integer,pointer:: items(:) => null() !serial numbers of lists items
    end type lsj_list

    ! Module state
    type(asf_basis_ls), public :: asf_set_ls
    integer, private :: limit_omitted_weight, ls_convention
    integer,pointer:: intermediate_two_j(:, :) ! xJ values
    real(real64) :: coefficient_threshold, max_omitted_percent

    ! Module constants
    character(len=1), private, parameter :: term_letters(0:20) =&
        [ "S","P","D","F","G","H","I","K","L","M","N","O","Q",   &
        "R","T","U","V","W","X","Y","Z" ]
    integer, parameter :: max_symmetry_blocks  = 20
    integer, parameter :: max_selected_levels = 100000
    character(*), private, parameter :: transformation_signature = "LSJTRN02", transformation_done = "LSJDONE1"

    ! Module state
    integer, private :: ls_csf_output_unit = -1
    integer, private :: lsj_label_unit = -1
    integer, private :: j_vector_output_unit = -1
    integer, private :: transformation_unit = -1
    integer, private :: unique_summary_unit = -1
    integer, private :: unique_label_unit = -1

contains
    ! Parallel rows are committed in the original LS order. The unchanged
    ! scalar backend remains the reference and handles sequential transformation-cache I/O.
    subroutine transform_asf_to_ls(parity_reference_position,retained_csf_mask,selected_level_count,symmetry_block_index, &
        & levels,first_csf_index, &
        & last_csf_index,transformation_output_mode)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real32
        use selected_level_indices, only: active_level_indices
!$    USE omp_lib, ONLY: omp_get_max_threads
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: parity_reference_position, selected_level_count, symmetry_block_index, first_csf_index, &
            & last_csf_index, transformation_output_mode
        integer, intent(in) :: retained_csf_mask(:), levels(max_symmetry_blocks,max_selected_levels)

        ! Local declarations
        integer :: threads, batch_size, first, last, row, ls, lev, level, status
        integer, allocatable :: selected(:)
        real(real64), allocatable :: rows(:,:), cumulative(:)
        logical, allocatable :: active(:)
        logical :: finished
        character(len=16) :: serial_backend

        threads = 1
!$    threads = omp_get_max_threads()
        call get_environment_variable("JJ2LSJ_SERIAL", serial_backend, status=status)
        if (status == 0) then
            if (trim(serial_backend) == "1") threads = 1
        end if
        if (threads <= 1 .or. transformation_output_mode /= 0) then
            call transform_asf_to_ls_serial(parity_reference_position,retained_csf_mask,selected_level_count, &
                & symmetry_block_index,levels,first_csf_index, &
                & last_csf_index,transformation_output_mode)
            return
        end if
        if (selected_level_count <= 0 .or. asf_set_ls%csf_set_ls%ls_csf_count <= 0) return

        ! Bound scratch to about 8 MiB, even for many selected levels. Each row
        ! is contiguous, with no per-pair allocations or cross-thread reduction.
        batch_size = min(256, 4*min(threads,64))
        batch_size = min(batch_size, max(1,1048576/selected_level_count))
        batch_size = min(batch_size, asf_set_ls%csf_set_ls%ls_csf_count)
        allocate(rows(selected_level_count,batch_size), active(batch_size))
        allocate(cumulative(selected_level_count), selected(selected_level_count))
        selected = levels(symmetry_block_index,1:selected_level_count)
        cumulative = 0.0_real64
        finished = .false.
        write(*,"(A,I0,A,I0)") " jj2lsj OpenMP: max threads=",threads, &
            ", LS rows per batch=",batch_size

        do first = 1, asf_set_ls%csf_set_ls%ls_csf_count, batch_size
            last = min(first+batch_size-1,asf_set_ls%csf_set_ls%ls_csf_count)
!$omp parallel do default(none) schedule(dynamic,1) &
!$omp shared(first,last,parity_reference_position,retained_csf_mask,selected,first_csf_index,last_csf_index,rows, &
!$omp& active) private(row,ls)
            do ls = first, last
                row = ls-first+1
                call transform_ls_row(ls,parity_reference_position,retained_csf_mask,selected,first_csf_index,last_csf_index, &
                    rows(:,row),active(row))
            end do
!$omp end parallel do

            ! Reproduce the scalar accumulation, strict comparison and first
            ! stopping row exactly; speculative later rows are never committed.
            do ls = first, last
                if (mod(ls,2000) == 0) print *, "LS_number=",ls
                row = ls-first+1
                if (.not. active(row)) cycle
                do lev = 1,selected_level_count
                    level = selected(lev)
                    cumulative(lev) = rows(lev,row)*rows(lev,row)+cumulative(lev)
                    asf_set_ls%asf(level)%eigenvector(ls) = rows(lev,row)
                    asf_set_ls%asf(level)%level_no = active_level_indices(level)
                end do
                if (limit_omitted_weight == 1) then
                    finished = .not. any(cumulative*100.0_real32+dabs(max_omitted_percent) <= 100.0_real32)
                    if (finished) exit
                end if
            end do
            if (finished) exit
        end do
        deallocate(rows,active,cumulative,selected)
    end subroutine transform_asf_to_ls

    recursive subroutine transform_ls_row(ls_csf_index,parity_reference_position,retained_csf_mask,selected, &
        first_csf_index,last_csf_index,values,active)
        ! Module dependencies
        use eigenstate_expansions, only: csf_mixing_coefficients
        use orbital_csf_basis, only: csf_count
        use csf_quantum_queries, only: csf_total_j_code, csf_parity_code
        use, intrinsic :: iso_fortran_env, only: int64
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: ls_csf_index, parity_reference_position, first_csf_index, last_csf_index
        integer, intent(in) :: retained_csf_mask(:), selected(:)
        real(real64), intent(out) :: values(:)
        logical, intent(out) :: active

        ! Local declarations
        integer :: jj_csf_index, selected_index, parity
        integer(int64) :: offset, mixing_column_stride
        real(real64) :: coefficient

        values = 0.0_real64
        parity = csf_parity_code(first_csf_index+parity_reference_position-1)
        active = (asf_set_ls%csf_set_ls%csf(ls_csf_index)%parity == "+" .and. parity == 1) .or. &
            (asf_set_ls%csf_set_ls%csf(ls_csf_index)%parity == "-" .and. parity == -1)
        if (.not. active) return
        mixing_column_stride = csf_count
        do jj_csf_index = first_csf_index, last_csf_index
            if (retained_csf_mask(jj_csf_index) /= 1) cycle
            if (asf_set_ls%csf_set_ls%csf(ls_csf_index)%two_j /= csf_total_j_code(jj_csf_index)-1) cycle
            coefficient = ls_jj_csf_overlap(jj_csf_index,ls_csf_index)
            do selected_index = 1,size(selected)
                offset = (selected(selected_index)-1)*mixing_column_stride
                values(selected_index) = values(selected_index)+csf_mixing_coefficients(jj_csf_index+offset)*coefficient
            end do
        end do
    end subroutine transform_ls_row

    ! Expands an atomic state functions from the same block,
    ! which is represented in a jj-coupling CSF basis into a basis
    ! of LS-coupling CSF.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: May 2011
    ! Modified by C. Cychen. 2021
    ! Modified by G. Gaigalas. 2022
    ! Modified by G. Gaigalas. 2024
    subroutine transform_asf_to_ls_serial(parity_reference_position,retained_csf_mask,selected_level_count, &
        & symmetry_block_index,levels,first_csf_index,     &
        last_csf_index,transformation_output_mode)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real32
        use, intrinsic :: iso_fortran_env, only: int64
        use eigenstate_expansions,          only: csf_mixing_coefficients
        use selected_level_indices,          only: active_level_indices
        use orbital_csf_basis,           only: csf_count
        use numeric_constants, only: zero
        use csf_quantum_queries, only: csf_total_j_code, csf_parity_code
        implicit none (type, external)

        ! Arguments and local declarations
        save  ! Preserve the legacy serial controller workspace.
        integer, intent(in) :: parity_reference_position, selected_level_count,symmetry_block_index,transformation_output_mode
        integer, intent(in) :: first_csf_index, last_csf_index
        integer, intent(in) :: retained_csf_mask(:)
        integer, intent(in) :: levels(max_symmetry_blocks, max_selected_levels)
        integer      :: jj_number, lev, level, ls_number, ijj, ils
        integer      :: weight_target_reached, cache_status
        real(real64) :: wa_transformation
        real(real64) :: ls_amplitudes(max_selected_levels)
        real(real64) :: cumulative_weights(max_selected_levels)
        integer(int64) :: level_mixing_offset, csf_count_int64
        cumulative_weights = zero
        csf_count_int64 = csf_count
        if (transformation_output_mode /= 0) call check_transformation_block(transformation_output_mode, &
            symmetry_block_index, retained_csf_mask(first_csf_index:last_csf_index))
        do ls_number = 1, asf_set_ls%csf_set_ls%ls_csf_count
            if(((ls_number/2000)*2000) == ls_number)                      &
                print*, "LS_number=",ls_number
            if ((asf_set_ls%csf_set_ls%csf(ls_number)%parity == "+"       &
                ! CYC-2023 Corrections, perform transfermation for both even and odd
                ! states
                ! .and.  ISPAR(iw1) == 1)  .or.                              &
                .and.  csf_parity_code(first_csf_index+parity_reference_position-1) == 1)  .or.                     &
                (asf_set_ls%csf_set_ls%csf(ls_number)%parity  == "-"       &
                ! .and.  ISPAR(iw1) == -1)) then
                .and.  csf_parity_code(first_csf_index+parity_reference_position-1) == -1)) then
                ! CYC-2023 end
                ls_amplitudes = zero
                do jj_number = first_csf_index, last_csf_index
                    if(retained_csf_mask(jj_number) == 1 .and.                        &
                        (asf_set_ls%csf_set_ls%csf(ls_number)%two_j ==       &
                        csf_total_j_code(jj_number)-1)) then
                        if(transformation_output_mode <= 1)                                       &
                            wa_transformation = ls_jj_csf_overlap(jj_number,ls_number)
                        if(transformation_output_mode == 1) write(transformation_unit)                             &
                            jj_number-first_csf_index+1,ls_number,wa_transformation
                        if(transformation_output_mode == 2) then
                            read(transformation_unit, iostat=cache_status) ijj, ils, wa_transformation
                            if (cache_status /= 0) error stop "Incomplete JJ2LSJ transformation cache"
                            if (ijj /= jj_number-first_csf_index+1 .or. ils /= ls_number) &
                                error stop "Unexpected coefficient indices in JJ2LSJ transformation cache"
                        end if
                        do lev = 1,selected_level_count
                            level = levels(symmetry_block_index,lev)
                            level_mixing_offset = (level-1)*csf_count_int64
                            ls_amplitudes(lev) =                                         &
                                ls_amplitudes(lev)+csf_mixing_coefficients(jj_number+level_mixing_offset)*wa_transformation
                        end do
                    end if
                end do
                weight_target_reached = 1
                do lev = 1,selected_level_count
                    level = levels(symmetry_block_index,lev)
                    cumulative_weights(lev) = ls_amplitudes(lev)*ls_amplitudes(lev) + cumulative_weights(lev)
                    asf_set_ls%asf(level)%eigenvector(ls_number) = ls_amplitudes(lev)
                    asf_set_ls%asf(level)%level_no = active_level_indices(level)
                    if((cumulative_weights(lev)*100.0_real32+dabs(max_omitted_percent)) <= 100.0_real32) weight_target_reached = 0
                end do
                if (limit_omitted_weight == 1 .and. weight_target_reached == 1) exit
            end if
        end do
        return
    end subroutine transform_asf_to_ls_serial

    subroutine check_transformation_block(mode, block, retained)
        integer, intent(in) :: mode, block, retained(:)
        integer :: status, stored_block
        integer, allocatable :: stored_mask(:)

        if (mode == 1) then
            write(transformation_unit) block, retained
        else
            allocate(stored_mask(size(retained)))
            read(transformation_unit, iostat=status) stored_block, stored_mask
            if (status /= 0) error stop "Incomplete JJ2LSJ transformation block"
            if (stored_block /= block .or. any(stored_mask /= retained)) &
                error stop "Selected blocks or retained CSFs differ from the JJ2LSJ transformation cache"
        end if
    end subroutine check_transformation_block

    !     Returns the value of the LS-jj transformation matrix for a given *
    !     set of quantum numbers form an array LS_jj_*_*.
    !
    !     Note that all (generalized) angular momentum quantum numbers
    !     except for l_shell are supplied as twice the physical value.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: May 2011
    recursive function ls_jj_recoupling_coefficient(l_shell, subshell_occupation, w, two_ls_quasispin, &
        & two_ls_term_l, two_ls_term_s, two_ls_term_j, lower_partner_orbital_two_j, lower_j_partner_occupation, &
        & lower_j_partner_two_quasispin, lower_j_partner_two_term_j, upper_partner_orbital_two_j, &
        & upper_j_partner_two_quasispin, upper_j_partner_two_term_j) result(coefficient_value)
        ! Module dependencies
        use numeric_constants,          only: zero, one
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: l_shell, subshell_occupation, w, two_ls_quasispin, two_ls_term_l, two_ls_term_s, two_ls_term_j
        integer, intent(in) :: lower_partner_orbital_two_j, lower_j_partner_occupation
        integer, intent(in) :: lower_j_partner_two_quasispin, lower_j_partner_two_term_j
        integer, intent(in) :: upper_partner_orbital_two_j, upper_j_partner_two_quasispin, upper_j_partner_two_term_j
        real(real64) :: coefficient_value

        ! Local declarations
        integer      :: nn, n1, n2, factor
        coefficient_value = zero; factor = one
        if (lower_partner_orbital_two_j < upper_partner_orbital_two_j .or. lower_partner_orbital_two_j == &
            & upper_partner_orbital_two_j)  then
            if (l_shell > 0 .and. subshell_occupation > 2*l_shell +1)  then
                nn = 4*l_shell + 2 - subshell_occupation;   n1 = lower_partner_orbital_two_j + 1 - lower_j_partner_occupation
                n2 = upper_partner_orbital_two_j + 1 - subshell_occupation + lower_j_partner_occupation
                if (mod(2*l_shell+1-two_ls_quasispin -((lower_partner_orbital_two_j+1)/2-lower_j_partner_two_quasispin)         &
                    -((upper_partner_orbital_two_j+1)/2-upper_j_partner_two_quasispin),4) /= 0) factor = - factor
            else
                nn = subshell_occupation;   n1 = lower_j_partner_occupation;   n2 = subshell_occupation - lower_j_partner_occupation
            end if
            if (nn == 1  .or.  nn == 0)  then
                if (nn == 0 .and. n1 == 0 .and. n2 == 0) then
                    coefficient_value = one
                else if (n1 == 1 .and. n2 == 0 .and. two_ls_term_j == lower_j_partner_two_term_j) then
                    coefficient_value = one
                else if (n1 == 0 .and. n2 == 1 .and. two_ls_term_j == upper_j_partner_two_term_j) then
                    coefficient_value = one
                else
                    coefficient_value = zero
                end if
            else if (nn == 2)  then
                if (two_ls_term_j > lower_j_partner_two_term_j + upper_j_partner_two_term_j  .or.  two_ls_term_j < &
                    & abs(lower_j_partner_two_term_j - upper_j_partner_two_term_j)) then
                    coefficient_value = zero
                else
                    if (n1 == 2 .and. n2 == 0) then
                        coefficient_value = single_partner_ls_jj_coefficient(l_shell,two_ls_term_l,two_ls_term_s,two_ls_term_j, &
                            & lower_partner_orbital_two_j,lower_partner_orbital_two_j)
                    else if (n1 == 1 .and. n2 == 1) then
                        coefficient_value = single_partner_ls_jj_coefficient(l_shell,two_ls_term_l,two_ls_term_s,two_ls_term_j, &
                            & lower_partner_orbital_two_j,upper_partner_orbital_two_j)
                    else if (n1 == 0 .and. n2 == 2) then
                        coefficient_value = single_partner_ls_jj_coefficient(l_shell,two_ls_term_l,two_ls_term_s,two_ls_term_j, &
                            & upper_partner_orbital_two_j,upper_partner_orbital_two_j)
                    end if
                end if
            else if (l_shell==1 .or. l_shell==2 .or. l_shell==3)  then
                coefficient_value = lookup_ls_jj_subshell_coefficient(l_shell,nn,w,two_ls_quasispin,two_ls_term_l,two_ls_term_s, &
                    & two_ls_term_j,n1,lower_j_partner_two_quasispin,lower_j_partner_two_term_j,upper_j_partner_two_quasispin, &
                        & upper_j_partner_two_term_j)
            end if
        else if (two_ls_term_l == 0)  then
            if (two_ls_term_s == two_ls_term_j .and. lower_partner_orbital_two_j == 1 .and. nn == n1)  then
                coefficient_value = one
            else if (two_ls_term_s == two_ls_term_j .and. upper_partner_orbital_two_j == 1 .and. nn == n2)  then
                coefficient_value = one
            else
                coefficient_value = zero
            end if
        else
            stop "ls_jj_recoupling_coefficient: program stop A."
        end if
        coefficient_value = factor * coefficient_value
    end function ls_jj_recoupling_coefficient

    !     Returns the value of the LS-jj transformation matrix
    !                                                (l^2 LSJ| j_1 j_2 J). *
    !
    !     Note that all (generalized) angular momentum quantum numbers
    !     except for l must be given as twice the original numbers, i.e.
    !     for the quantum numbers L, S, J, jm_shell, jp_shell.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: May 2011
    recursive function single_partner_ls_jj_coefficient(l_shell,l,s,j,jm_shell,jp_shell)   result(coefficient_value)
        ! Module dependencies
        use numeric_constants,          only: zero, one
        use wigner_9j, only: wigner_9j_value
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: l_shell, l, s, j, jm_shell, jp_shell
        real(real64)        :: coefficient_value

        ! Local declarations
        integer      :: delta_j
        real(real64) :: rac9
        coefficient_value = zero
        if (mod(l+s,4) /= 0)  then
            coefficient_value = zero
        else if (mod(l+s+j,2) /= 0)  then
            coefficient_value = zero
        else if (mod(jm_shell+jp_shell+j,2) /= 0)  then
            coefficient_value = zero
        else
            if (jm_shell == jp_shell) then
                if (mod(j,4) /= 0)  then
                    coefficient_value = zero
                else
                    coefficient_value = (jm_shell+one) * sqrt((l+one) * (s+one))
                end if
            else
                coefficient_value = sqrt(2*(l+one)*(s+one)*(jm_shell+one)*(jp_shell+one))
            end if
            call wigner_9j_value(2*l_shell,2*l_shell,l,1,1,s,jm_shell,jp_shell,j,   &
                1,delta_j,rac9)
            if (delta_j /= 0)  then
                call wigner_9j_value(2*l_shell,2*l_shell,l,1,1,s,jm_shell,jp_shell,j,&
                    0,delta_j,rac9)
                coefficient_value = coefficient_value * rac9
            end if
        end if
    end function single_partner_ls_jj_coefficient

    !     Returns the value of the LS-jj transformation matrix for a given *
    !     set of quantum numbers.
    !
    !     Note that all (generalized) angular momentum quantum numbers
    !     except for l must be given as twice the original numbers, i.e.
    !     for the quantum numbers Q, L, S, J, Qm, Jm, Qp, Jp.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: May 2011
    recursive function lookup_ls_jj_subshell_coefficient(lshell,subshell_occupation,w,two_ls_quasispin,two_ls_term_l, &
        & two_ls_term_s,two_ls_term_j,lower_j_partner_occupation,lower_j_partner_two_quasispin,lower_j_partner_two_term_j, &
            & upper_j_partner_two_quasispin,upper_j_partner_two_term_j)  result(coefficient_value)
        ! Module dependencies
        use ls_jj_subshell_terms, only: ls_jj_number_d10, ls_jj_number_d3, ls_jj_number_d4, ls_jj_number_d5, ls_jj_number_d6, &
            & ls_jj_number_d7, ls_jj_number_d8, ls_jj_number_d9, ls_jj_number_f3, ls_jj_number_f4, ls_jj_number_f5, &
            & ls_jj_number_f6, ls_jj_number_f7, ls_jj_number_p3, ls_jj_number_p4, ls_jj_number_p5, ls_jj_number_p6
        use ls_jj_coefficients_low_j, only: ls_jj_d_10, ls_jj_d_3, ls_jj_d_4, ls_jj_d_5, ls_jj_d_6, ls_jj_d_7, ls_jj_d_8, &
            & ls_jj_d_9, ls_jj_f_3, ls_jj_f_4, ls_jj_f_5, ls_jj_number_d10, ls_jj_number_d3, ls_jj_number_d4, &
            & ls_jj_number_d5, ls_jj_number_d6, ls_jj_number_d7, ls_jj_number_d8, ls_jj_number_d9, ls_jj_number_f3, &
            & ls_jj_number_f4, ls_jj_number_f5, ls_jj_number_f6, ls_jj_number_f7, ls_jj_number_p3, ls_jj_number_p4, &
            & ls_jj_number_p5, ls_jj_number_p6, ls_jj_p_3, ls_jj_p_4, ls_jj_p_5, ls_jj_p_6
        use ls_jj_coefficients_f6, only: ls_jj_f_6, ls_jj_number_d10, ls_jj_number_d3, ls_jj_number_d4, ls_jj_number_d5, &
            & ls_jj_number_d6, ls_jj_number_d7, ls_jj_number_d8, ls_jj_number_d9, ls_jj_number_f3, ls_jj_number_f4, &
            & ls_jj_number_f5, ls_jj_number_f6, ls_jj_number_f7, ls_jj_number_p3, ls_jj_number_p4, ls_jj_number_p5, &
            & ls_jj_number_p6
        use ls_jj_coefficients_f7, only: ls_jj_f_7, ls_jj_number_d10, ls_jj_number_d3, ls_jj_number_d4, ls_jj_number_d5, &
            & ls_jj_number_d6, ls_jj_number_d7, ls_jj_number_d8, ls_jj_number_d9, ls_jj_number_f3, ls_jj_number_f4, &
            & ls_jj_number_f5, ls_jj_number_f6, ls_jj_number_f7, ls_jj_number_p3, ls_jj_number_p4, ls_jj_number_p5, &
            & ls_jj_number_p6
        use numeric_constants,          only: zero, one
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: lshell, subshell_occupation, w, two_ls_quasispin, two_ls_term_l, two_ls_term_s, two_ls_term_j, &
            & lower_j_partner_occupation, lower_j_partner_two_quasispin, lower_j_partner_two_term_j, &
                & upper_j_partner_two_quasispin, upper_j_partner_two_term_j
        real(real64)        :: coefficient_value

        ! Local declarations
        integer :: i
        coefficient_value = zero
        if (lshell == 0) then
        else if (lshell == 1) then
            select case(subshell_occupation)
            case(3)
                ! Use data from the array LS_jj_p_3
                do  i = 1,ls_jj_number_p3
                    if(w  ==ls_jj_p_3(i)%w  .and. two_ls_quasispin  ==ls_jj_p_3(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_p_3(i)%l  .and. two_ls_term_s  ==ls_jj_p_3(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_p_3(i)%j  .and. lower_j_partner_occupation ==ls_jj_p_3(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_p_3(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_p_3(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_p_3(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_p_3(i)%jp) then
                        coefficient_value = ls_jj_p_3(i)%factor * &
                            sqrt( one*ls_jj_p_3(i)%nom/ls_jj_p_3(i)%denom )
                        return
                    end if
                end do
            case(4)
                ! Use data from the array LS_jj_p_4
                do  i = 1,ls_jj_number_p4
                    if(w  ==ls_jj_p_4(i)%w  .and. two_ls_quasispin  ==ls_jj_p_4(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_p_4(i)%l  .and. two_ls_term_s  ==ls_jj_p_4(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_p_4(i)%j  .and. lower_j_partner_occupation ==ls_jj_p_4(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_p_4(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_p_4(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_p_4(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_p_4(i)%jp) then
                        coefficient_value = ls_jj_p_4(i)%factor * &
                            sqrt( one*ls_jj_p_4(i)%nom/ls_jj_p_4(i)%denom )
                        return
                    end if
                end do
            case(5)
                ! Use data from the array LS_jj_p_5
                do  i = 1,ls_jj_number_p5
                    if(w  ==ls_jj_p_5(i)%w  .and. two_ls_quasispin  ==ls_jj_p_5(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_p_5(i)%l  .and. two_ls_term_s  ==ls_jj_p_5(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_p_5(i)%j  .and. lower_j_partner_occupation ==ls_jj_p_5(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_p_5(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_p_5(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_p_5(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_p_5(i)%jp) then
                        coefficient_value = ls_jj_p_5(i)%factor * &
                            sqrt( one*ls_jj_p_5(i)%nom/ls_jj_p_5(i)%denom )
                        return
                    end if
                end do
            case(6)
                ! Use data from the array LS_jj_p_6
                do  i = 1,ls_jj_number_p6
                    if(w  ==ls_jj_p_6(i)%w  .and. two_ls_quasispin  ==ls_jj_p_6(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_p_6(i)%l  .and. two_ls_term_s  ==ls_jj_p_6(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_p_6(i)%j  .and. lower_j_partner_occupation ==ls_jj_p_6(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_p_6(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_p_6(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_p_6(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_p_6(i)%jp) then
                        coefficient_value = ls_jj_p_6(i)%factor * &
                            sqrt( one*ls_jj_p_6(i)%nom/ls_jj_p_6(i)%denom )
                        return
                    end if
                end do
            case default
                stop "lookup_ls_jj_subshell_coefficient: program stop A."
            end select
        else if (lshell == 2) then
            select case(subshell_occupation)
            case(3)
                ! Use data from the array LS_jj_d_3
                do  i = 1,ls_jj_number_d3
                    if(w  ==ls_jj_d_3(i)%w  .and. two_ls_quasispin  ==ls_jj_d_3(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_d_3(i)%l  .and. two_ls_term_s  ==ls_jj_d_3(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_d_3(i)%j  .and. lower_j_partner_occupation ==ls_jj_d_3(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_d_3(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_d_3(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_d_3(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_d_3(i)%jp) then
                        coefficient_value = ls_jj_d_3(i)%factor * &
                            sqrt( one*ls_jj_d_3(i)%nom/ls_jj_d_3(i)%denom )
                        return
                    end if
                end do
            case(4)
                ! Use data from the array LS_jj_d_4
                do  i = 1,ls_jj_number_d4
                    if(w  ==ls_jj_d_4(i)%w  .and. two_ls_quasispin  ==ls_jj_d_4(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_d_4(i)%l  .and. two_ls_term_s  ==ls_jj_d_4(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_d_4(i)%j  .and. lower_j_partner_occupation ==ls_jj_d_4(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_d_4(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_d_4(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_d_4(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_d_4(i)%jp) then
                        coefficient_value = ls_jj_d_4(i)%factor * &
                            sqrt( one*ls_jj_d_4(i)%nom/ls_jj_d_4(i)%denom )
                        return
                    end if
                end do
            case(5)
                ! Use data from the array LS_jj_d_5
                do  i = 1,ls_jj_number_d5
                    if(w  ==ls_jj_d_5(i)%w  .and. two_ls_quasispin  ==ls_jj_d_5(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_d_5(i)%l  .and. two_ls_term_s  ==ls_jj_d_5(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_d_5(i)%j  .and. lower_j_partner_occupation ==ls_jj_d_5(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_d_5(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_d_5(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_d_5(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_d_5(i)%jp) then
                        coefficient_value = ls_jj_d_5(i)%factor * &
                            sqrt( one*ls_jj_d_5(i)%nom/ls_jj_d_5(i)%denom )
                        return
                    end if
                end do
            case(6)
                ! Use data from the array LS_jj_d_6
                do  i = 1,ls_jj_number_d6
                    if(w  ==ls_jj_d_6(i)%w  .and. two_ls_quasispin  ==ls_jj_d_6(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_d_6(i)%l  .and. two_ls_term_s  ==ls_jj_d_6(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_d_6(i)%j  .and. lower_j_partner_occupation ==ls_jj_d_6(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_d_6(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_d_6(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_d_6(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_d_6(i)%jp) then
                        coefficient_value = ls_jj_d_6(i)%factor * &
                            sqrt( one*ls_jj_d_6(i)%nom/ls_jj_d_6(i)%denom )
                        return
                    end if
                end do
            case(7)
                ! Use data from the array LS_jj_d_7
                do  i = 1,ls_jj_number_d7
                    if(w  ==ls_jj_d_7(i)%w  .and. two_ls_quasispin  ==ls_jj_d_7(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_d_7(i)%l  .and. two_ls_term_s  ==ls_jj_d_7(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_d_7(i)%j  .and. lower_j_partner_occupation ==ls_jj_d_7(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_d_7(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_d_7(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_d_7(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_d_7(i)%jp) then
                        coefficient_value = ls_jj_d_7(i)%factor * &
                            sqrt( one*ls_jj_d_7(i)%nom/ls_jj_d_7(i)%denom )
                        return
                    end if
                end do
            case(8)
                ! Use data from the array LS_jj_d_8
                do  i = 1,ls_jj_number_d8
                    if(w  ==ls_jj_d_8(i)%w  .and. two_ls_quasispin  ==ls_jj_d_8(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_d_8(i)%l  .and. two_ls_term_s  ==ls_jj_d_8(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_d_8(i)%j  .and. lower_j_partner_occupation ==ls_jj_d_8(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_d_8(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_d_8(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_d_8(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_d_8(i)%jp) then
                        coefficient_value = ls_jj_d_8(i)%factor * &
                            sqrt( one*ls_jj_d_8(i)%nom/ls_jj_d_8(i)%denom )
                        return
                    end if
                end do
            case(9)
                ! Use data from the array LS_jj_d_9
                do  i = 1,ls_jj_number_d9
                    if(w  ==ls_jj_d_9(i)%w  .and. two_ls_quasispin  ==ls_jj_d_9(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_d_9(i)%l  .and. two_ls_term_s  ==ls_jj_d_9(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_d_9(i)%j  .and. lower_j_partner_occupation ==ls_jj_d_9(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_d_9(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_d_9(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_d_9(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_d_9(i)%jp) then
                        coefficient_value = ls_jj_d_9(i)%factor * &
                            sqrt( one*ls_jj_d_9(i)%nom/ls_jj_d_9(i)%denom )
                        return
                    end if
                end do
            case(10)
                ! Use data from the array LS_jj_d_10
                do  i = 1,ls_jj_number_d10
                    if(w ==ls_jj_d_10(i)%w  .and. two_ls_quasispin ==ls_jj_d_10(i)%q  .and. &
                        two_ls_term_l ==ls_jj_d_10(i)%l  .and. two_ls_term_s ==ls_jj_d_10(i)%s  .and. &
                        two_ls_term_j ==ls_jj_d_10(i)%j  .and. lower_j_partner_occupation==ls_jj_d_10(i)%nm .and. &
                        lower_j_partner_two_quasispin==ls_jj_d_10(i)%qm .and. lower_j_partner_two_term_j==ls_jj_d_10(i)%jm .and. &
                        upper_j_partner_two_quasispin==ls_jj_d_10(i)%qp .and. upper_j_partner_two_term_j==ls_jj_d_10(i)%jp) then
                        coefficient_value= ls_jj_d_10(i)%factor * &
                            sqrt( one*ls_jj_d_10(i)%nom/ls_jj_d_10(i)%denom )
                        return
                    end if
                end do
            case default
                stop "lookup_ls_jj_subshell_coefficient: program stop B."
            end select
        else if (lshell == 3) then
            select case(subshell_occupation)
            case(3)
                ! Use data from the array LS_jj_f_3
                do  i = 1,ls_jj_number_f3
                    if(w  ==ls_jj_f_3(i)%w  .and. two_ls_quasispin  ==ls_jj_f_3(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_f_3(i)%l  .and. two_ls_term_s  ==ls_jj_f_3(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_f_3(i)%j  .and. lower_j_partner_occupation ==ls_jj_f_3(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_f_3(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_f_3(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_f_3(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_f_3(i)%jp) then
                        coefficient_value = ls_jj_f_3(i)%factor * &
                            sqrt( one*ls_jj_f_3(i)%nom/ls_jj_f_3(i)%denom )
                        return
                    end if
                end do
            case(4)
                ! Use data from the array LS_jj_f_4
                do  i = 1,ls_jj_number_f4
                    if(w  ==ls_jj_f_4(i)%w  .and. two_ls_quasispin  ==ls_jj_f_4(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_f_4(i)%l  .and. two_ls_term_s  ==ls_jj_f_4(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_f_4(i)%j  .and. lower_j_partner_occupation ==ls_jj_f_4(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_f_4(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_f_4(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_f_4(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_f_4(i)%jp) then
                        coefficient_value = ls_jj_f_4(i)%factor * &
                            sqrt( one*ls_jj_f_4(i)%nom/ls_jj_f_4(i)%denom )
                        return
                    end if
                end do
            case(5)
                ! Use data from the array LS_jj_f_5
                do  i = 1,ls_jj_number_f5
                    if(w  ==ls_jj_f_5(i)%w  .and. two_ls_quasispin  ==ls_jj_f_5(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_f_5(i)%l  .and. two_ls_term_s  ==ls_jj_f_5(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_f_5(i)%j  .and. lower_j_partner_occupation ==ls_jj_f_5(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_f_5(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_f_5(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_f_5(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_f_5(i)%jp) then
                        coefficient_value = ls_jj_f_5(i)%factor * &
                            sqrt( one*ls_jj_f_5(i)%nom/ls_jj_f_5(i)%denom )
                        return
                    end if
                end do
            case(6)
                ! Use data from the array LS_jj_f_6
                do  i = 1,ls_jj_number_f6
                    if(w  ==ls_jj_f_6(i)%w  .and. two_ls_quasispin  ==ls_jj_f_6(i)%q  .and. &
                        two_ls_term_l  ==ls_jj_f_6(i)%l  .and. two_ls_term_s  ==ls_jj_f_6(i)%s  .and. &
                        two_ls_term_j  ==ls_jj_f_6(i)%j  .and. lower_j_partner_occupation ==ls_jj_f_6(i)%nm .and. &
                        lower_j_partner_two_quasispin ==ls_jj_f_6(i)%qm .and. lower_j_partner_two_term_j ==ls_jj_f_6(i)%jm .and. &
                        upper_j_partner_two_quasispin ==ls_jj_f_6(i)%qp .and. upper_j_partner_two_term_j ==ls_jj_f_6(i)%jp) then
                        coefficient_value = ls_jj_f_6(i)%factor * &
                            sqrt( one*ls_jj_f_6(i)%nom/ls_jj_f_6(i)%denom )
                        return
                    end if
                end do
            case(7)
                ! Use data from the array LS_jj_f_7
                do  i = 1, ls_jj_number_f7, 1
                    if (w == ls_jj_f_7(i)%w) then
                        if(two_ls_quasispin == ls_jj_f_7(i)%q) then
                            if(two_ls_term_l == ls_jj_f_7(i)%l) then
                                if(two_ls_term_s == ls_jj_f_7(i)%s) then
                                    if(two_ls_term_j == ls_jj_f_7(i)%j) then
                                        if(lower_j_partner_occupation == ls_jj_f_7(i)%nm) then
                                            if(lower_j_partner_two_quasispin == ls_jj_f_7(i)%qm .and. &
                                                lower_j_partner_two_term_j == ls_jj_f_7(i)%jm .and. &
                                                upper_j_partner_two_quasispin == ls_jj_f_7(i)%qp .and. &
                                                upper_j_partner_two_term_j == ls_jj_f_7(i)%jp) then
                                                coefficient_value = one * ls_jj_f_7(i)%factor  * &
                                                    dsqrt(one*ls_jj_f_7(i)%nom)/ &
                                                    dsqrt(one*ls_jj_f_7(i)%denom )
                                                return
                                            end if
                                        end if
                                    end if
                                end if
                            end if
                        end if
                    end if
                end do
            case default
                stop "lookup_ls_jj_subshell_coefficient: program stop C."
            end select
        else
            stop "lookup_ls_jj_subshell_coefficient: program stop D."
        end if
    end function lookup_ls_jj_subshell_coefficient

    !     Dellocates the storage of asf_set_LS.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: May 2011
    subroutine release_ls_expansion(asf_set_ls)
        ! Module dependencies
        use selected_level_indices,             only: active_level_count
        use csf_quantum_queries, only: csf_total_j_code
        implicit none (type, external)

        ! Arguments and local declarations
        save  ! Preserve the legacy serial controller workspace.
        type(asf_basis_ls), intent(inout) :: asf_set_ls
        integer :: i
        deallocate(asf_set_ls%csf_set_ls%shell)
        deallocate(asf_set_ls%csf_set_ls%parent)
        do i = 1, asf_set_ls%csf_set_ls%ls_csf_count, 1
            deallocate(asf_set_ls%csf_set_ls%csf(i)%occupation)
            deallocate(asf_set_ls%csf_set_ls%csf(i)%seniority)
            deallocate(asf_set_ls%csf_set_ls%csf(i)%w)
            deallocate(asf_set_ls%csf_set_ls%csf(i)%shell_two_l)
            deallocate(asf_set_ls%csf_set_ls%csf(i)%shell_two_s)
            deallocate(asf_set_ls%csf_set_ls%csf(i)%intermediate_two_l)
            deallocate(asf_set_ls%csf_set_ls%csf(i)%intermediate_two_s)
        end do
        deallocate(asf_set_ls%csf_set_ls%csf)
        do i = 1, active_level_count
            deallocate(asf_set_ls%asf(i)%eigenvector)
        end do
        deallocate(asf_set_ls%asf)
    end subroutine release_ls_expansion

    !     The spectroscopic notation of a shell in LS coupling is returned.*
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: May 2011
    subroutine format_ls_csf(csf_index,shell_number,string_ls, string_xls)
        ! Module dependencies
        use half_integer_format, only: format_half_integer
        implicit none (type, external)

        ! Arguments and local declarations
        save  ! Preserve the legacy serial controller workspace.
        integer, intent(in)           :: csf_index,shell_number
        character(len=4), intent(out) :: string_ls, string_xls
        integer          :: text_length
        character(len=1) :: string_s, string_v
        character(len=4) :: string_cnum
        call                                                            &
            format_half_integer(2*(1+asf_set_ls%csf_set_ls%csf(csf_index)%shell_two_s(shell_number)),&
            string_cnum,text_length)
        string_s = string_cnum(1:text_length)
        if (asf_set_ls%csf_set_ls%shell(shell_number)%l < 3)    then
            call                                                         &
                format_half_integer(2*asf_set_ls%csf_set_ls%csf(csf_index)%seniority(shell_number),&
                string_cnum,text_length)
            string_v = string_cnum(1:text_length)
        else
            call                                                         &
                format_half_integer(2*asf_set_ls%csf_set_ls%csf(csf_index)%w(shell_number),&
                string_cnum,text_length)
            string_v = string_cnum(1:text_length)
        end if
        string_ls = string_s //                                         &
            term_letters(asf_set_ls%csf_set_ls%csf(csf_index)%shell_two_l(shell_number)/2)&
            //string_v
        call                                                            &
            format_half_integer(2*(1+asf_set_ls%csf_set_ls%csf(csf_index)%intermediate_two_s(shell_number)),&
            string_cnum,text_length)
        string_s =  string_cnum(1:text_length)
        string_xls = string_s //                                        &
            term_letters(asf_set_ls%csf_set_ls%csf(csf_index)%intermediate_two_l(shell_number)/2)
    end subroutine format_ls_csf

    !     The spectroscopic notation of a shell in LS coupling is returned.*
    !
    ! Written by G. Gaigalas,
    ! Modified by G. Gaigalas. 2024
    subroutine format_ls_csf_label(csf_index,shell_number,string_ls,string_xls)
        ! Module dependencies
        use half_integer_format, only: format_half_integer
        implicit none (type, external)

        ! Arguments and local declarations
        save  ! Preserve the legacy serial controller workspace.
        integer, intent(in)           :: csf_index,shell_number
        character(len=4), intent(out) :: string_ls, string_xls
        integer          :: text_length, w, tl, ts, seniority
        character(len=1) :: string_s, string_v
        character(len=4) :: string_cnum
        call                                                            &
            format_half_integer(2*(1+asf_set_ls%csf_set_ls%csf(csf_index)%shell_two_s(shell_number)),&
            string_cnum,text_length)
        string_s = string_cnum(1:text_length)
        ! 2024
        if (asf_set_ls%csf_set_ls%shell(shell_number)%l == 2)    then
            if (ls_convention == 0) then
                call                                                       &
                    format_half_integer(2*asf_set_ls%csf_set_ls%csf(csf_index)% &
                    seniority(shell_number),string_cnum,text_length)
                string_v = string_cnum(1:text_length)
            else
                ts=asf_set_ls%csf_set_ls%csf(csf_index)%shell_two_s(shell_number)
                tl=asf_set_ls%csf_set_ls%csf(csf_index)%shell_two_l(shell_number)
                seniority=asf_set_ls%csf_set_ls%csf(csf_index)%seniority(shell_number)
                if(ts+1 == 1) then
                    if(tl == 0) then
                        ! 1S
                        if(seniority == 0) then
                            w = 1
                        else if(seniority == 4) then
                            w = 2
                        else
                            stop "format_ls_csf_label(): program stop 1S."
                        end if
                    else if(tl == 4) then
                        ! 1D
                        if(seniority == 2) then
                            w = 1
                        else if(seniority == 4) then
                            w = 2
                        else
                            stop "format_ls_csf_label(): program stop 1D."
                        end if
                    else if(tl == 8) then
                        ! 1G
                        if(seniority == 2) then
                            w = 1
                        else if(seniority == 4) then
                            w = 2
                        else
                            stop "format_ls_csf_label(): program stop 1G."
                        end if
                    else
                        w = 0
                    end if
                else if(ts+1 == 2) then
                    if(tl == 4) then
                        ! 2D
                        if(seniority == 1) then
                            w = 1
                        else if(seniority == 3) then
                            w = 2
                        else if(seniority == 5) then
                            w = 3
                        else
                            stop "format_ls_csf_label(): program stop 2D."
                        end if
                    else if(tl == 6) then
                        ! 2F
                        if(seniority == 3) then
                            w = 1
                        else if(seniority == 5) then
                            w = 2
                        else
                            stop "format_ls_csf_label(): program stop 2F."
                        end if
                    else if(tl == 8) then
                        ! 2G
                        if(seniority == 3) then
                            w = 1
                        else if(seniority == 5) then
                            w = 2
                        else
                            stop "format_ls_csf_label(): program stop 2G."
                        end if
                    else
                        w = 0
                    end if
                else if(ts+1 == 3) then
                    if(tl == 2) then
                        ! 3P
                        if(seniority == 2) then
                            w = 1
                        else if(seniority == 4) then
                            w = 2
                        else
                            stop "format_ls_csf_label(): program stop 3P."
                        end if
                    else if(tl == 6) then
                        ! 3F
                        if(seniority == 2) then
                            w = 1
                        else if(seniority == 4) then
                            w = 2
                        else
                            stop "format_ls_csf_label(): program stop 3F."
                        end if
                    else
                        w = 0
                    end if
                else
                    w = 0
                end if
                call format_half_integer(2*w,string_cnum,text_length)
                string_v = string_cnum(1:text_length)
            end if
        else if (asf_set_ls%csf_set_ls%shell(shell_number)%l < 3)    then
            if (ls_convention == 0) then
                call                                                       &
                    format_half_integer(2*asf_set_ls%csf_set_ls%csf(csf_index)% &
                    seniority(shell_number),string_cnum,text_length)
                string_v = string_cnum(1:text_length)
            else
                w = asf_set_ls%csf_set_ls%csf(csf_index)%w(shell_number)
                call format_half_integer(2*w,string_cnum,text_length)
                string_v = string_cnum(1:text_length)
            end if
        else
            w = asf_set_ls%csf_set_ls%csf(csf_index)%w(shell_number)
            call format_half_integer(2*w, string_cnum,text_length)
            string_v = string_cnum(1:text_length)
        end if
        if(asf_set_ls%csf_set_ls%shell(shell_number)%l < 2) then
            string_ls = string_s //                                      &
                term_letters(asf_set_ls%csf_set_ls%csf(csf_index)%              &
                shell_two_l(shell_number)/2)
        else if(asf_set_ls%csf_set_ls%shell(shell_number)%l == 3        &
            .and. w == 0) then
            string_ls = string_s //                                      &
                term_letters(asf_set_ls%csf_set_ls%csf(csf_index)%              &
                shell_two_l(shell_number)/2)
        else if(ls_convention == 1 .and. w == 0) then
            string_ls = string_s //                                      &
                term_letters(asf_set_ls%csf_set_ls%csf(csf_index)%              &
                shell_two_l(shell_number)/2)
        else
            if(w == 10) then
                string_ls = string_s //                                   &
                    term_letters(asf_set_ls%csf_set_ls%csf(csf_index)%shell_two_l(    &
                    shell_number)/2)//"A"
            else
                string_ls = string_s //                                   &
                    term_letters(asf_set_ls%csf_set_ls%csf(csf_index)%shell_two_l(    &
                    shell_number)/2)//string_v
            end if
        end if
        ! 2024 end
        call                                                            &
            format_half_integer(2*(1+asf_set_ls%csf_set_ls%csf(csf_index)%intermediate_two_s(shell_number)),&
            string_cnum,text_length)
        string_s =  string_cnum(1:text_length)
        string_xls = string_s //                                        &
            term_letters(asf_set_ls%csf_set_ls%csf(csf_index)%intermediate_two_l(shell_number)/2)
    end subroutine format_ls_csf_label

    subroutine select_ls_transition_labels
        use orbital_csf_basis, only: csf_count, orbital_count
        use csf_quantum_queries, only: csf_cumulative_coupling
        implicit none (type, external)
        integer :: csf_index, orbital_index

        do csf_index = 1, csf_count
            do orbital_index = 1, orbital_count
                intermediate_two_j(orbital_index,csf_index) = csf_cumulative_coupling(orbital_index,csf_index) - 1
            end do
        end do
    end subroutine  select_ls_transition_labels

    !     This procedure returns all allowed subshell terms (l,w,Q,L,S)
    !     for given l^N, N = 0, 1, 2 or 3.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: May 2011
    subroutine enumerate_ls_subshell_terms(l_shell,n,ls,number)
        ! Module dependencies
        use ls_jj_subshell_terms, only: subshell_term_ls, term_ls_d, term_ls_f, term_ls_g1, term_ls_g2, term_ls_h1, term_ls_h2, &
            & term_ls_i1, term_ls_i2, term_ls_k1, term_ls_k2, term_ls_l1, term_ls_l2, &
            & term_ls_p, term_ls_s
        implicit none (type, external)

        ! Arguments and local declarations
        save  ! Preserve the legacy serial controller workspace.
        integer, intent(in)                                 :: l_shell, n
        type(subshell_term_ls), intent(out) :: ls(120)
        integer, intent(out)                                :: number
        integer :: m_q, i, j
        m_q = n - 2* l_shell - 1;  j = 0
        select case (l_shell)
        case (0)
            do i = 1,2
                if (mod(m_q + term_ls_s(i)%subshell_two_quasispin,2) == 0) then
                    if (abs(m_q) <= term_ls_s(i)%subshell_two_quasispin) then
                        j = j + 1
                        ls(j)%subshell_l  = term_ls_s(i)%subshell_l
                        ls(j)%w        = term_ls_s(i)%w
                        ls(j)%subshell_two_quasispin        = term_ls_s(i)%subshell_two_quasispin
                        ls(j)%subshell_two_l       = term_ls_s(i)%subshell_two_l
                        ls(j)%subshell_two_s        = term_ls_s(i)%subshell_two_s
                    end if
                end if
            end do
        case (1)
            do i = 1,6
                if (mod(m_q + term_ls_p(i)%subshell_two_quasispin,2) == 0) then
                    if (abs(m_q) <= term_ls_p(i)%subshell_two_quasispin) then
                        j = j + 1
                        ls(j)%subshell_l  = term_ls_p(i)%subshell_l
                        ls(j)%w        = term_ls_p(i)%w
                        ls(j)%subshell_two_quasispin        = term_ls_p(i)%subshell_two_quasispin
                        ls(j)%subshell_two_l       = term_ls_p(i)%subshell_two_l
                        ls(j)%subshell_two_s        = term_ls_p(i)%subshell_two_s
                    end if
                end if
            end do
        case (2)
            do i = 1,32
                if (mod(m_q + term_ls_d(i)%subshell_two_quasispin,2) == 0) then
                    if (abs(m_q) <= term_ls_d(i)%subshell_two_quasispin) then
                        j = j + 1
                        ls(j)%subshell_l  = term_ls_d(i)%subshell_l
                        ls(j)%w        = term_ls_d(i)%w
                        ls(j)%subshell_two_quasispin        = term_ls_d(i)%subshell_two_quasispin
                        ls(j)%subshell_two_l       = term_ls_d(i)%subshell_two_l
                        ls(j)%subshell_two_s        = term_ls_d(i)%subshell_two_s
                    end if
                end if
            end do
        case (3)
            do i = 1,238
                if (mod(m_q + term_ls_f(i)%subshell_two_quasispin,2) == 0) then
                    if (abs(m_q) <= term_ls_f(i)%subshell_two_quasispin) then
                        j = j + 1
                        ls(j)%subshell_l  = term_ls_f(i)%subshell_l
                        ls(j)%w        = term_ls_f(i)%w
                        ls(j)%subshell_two_quasispin        = term_ls_f(i)%subshell_two_quasispin
                        ls(j)%subshell_two_l       = term_ls_f(i)%subshell_two_l
                        ls(j)%subshell_two_s        = term_ls_f(i)%subshell_two_s
                    end if
                end if
            end do
        case (4)
            select case (n)
            case (1)
                i = 1;       j = 1
                ls(j)%subshell_l  = term_ls_g1(i)%subshell_l
                ls(j)%w        = term_ls_g1(i)%w
                ls(j)%subshell_two_quasispin        = term_ls_g1(i)%subshell_two_quasispin
                ls(j)%subshell_two_l       = term_ls_g1(i)%subshell_two_l
                ls(j)%subshell_two_s        = term_ls_g1(i)%subshell_two_s
            case (2)
                do i = 1,9
                    j = j + 1
                    ls(j)%subshell_l  = term_ls_g2(i)%subshell_l
                    ls(j)%w        = term_ls_g2(i)%w
                    ls(j)%subshell_two_quasispin        = term_ls_g2(i)%subshell_two_quasispin
                    ls(j)%subshell_two_l       = term_ls_g2(i)%subshell_two_l
                    ls(j)%subshell_two_s        = term_ls_g2(i)%subshell_two_s
                end do
            case default
                stop  "enumerate_ls_subshell_terms(): program stop A."
            end select
        case (5)
            select case (n)
            case (1)
                i = 1;       j = 1
                ls(j)%subshell_l  = term_ls_h1(i)%subshell_l
                ls(j)%w        = term_ls_h1(i)%w
                ls(j)%subshell_two_quasispin        = term_ls_h1(i)%subshell_two_quasispin
                ls(j)%subshell_two_l       = term_ls_h1(i)%subshell_two_l
                ls(j)%subshell_two_s        = term_ls_h1(i)%subshell_two_s
            case (2)
                do i = 1,11
                    j = j + 1
                    ls(j)%subshell_l  = term_ls_h2(i)%subshell_l
                    ls(j)%w        = term_ls_h2(i)%w
                    ls(j)%subshell_two_quasispin        = term_ls_h2(i)%subshell_two_quasispin
                    ls(j)%subshell_two_l       = term_ls_h2(i)%subshell_two_l
                    ls(j)%subshell_two_s        = term_ls_h2(i)%subshell_two_s
                end do
            case default
                stop  "enumerate_ls_subshell_terms(): program stop B."
            end select
        case (6)
            select case (n)
            case (1)
                i = 1;       j = 1
                ls(j)%subshell_l  = term_ls_i1(i)%subshell_l
                ls(j)%w        = term_ls_i1(i)%w
                ls(j)%subshell_two_quasispin        = term_ls_i1(i)%subshell_two_quasispin
                ls(j)%subshell_two_l       = term_ls_i1(i)%subshell_two_l
                ls(j)%subshell_two_s        = term_ls_i1(i)%subshell_two_s
            case (2)
                do i = 1,13
                    j = j + 1
                    ls(j)%subshell_l  = term_ls_i2(i)%subshell_l
                    ls(j)%w        = term_ls_i2(i)%w
                    ls(j)%subshell_two_quasispin        = term_ls_i2(i)%subshell_two_quasispin
                    ls(j)%subshell_two_l       = term_ls_i2(i)%subshell_two_l
                    ls(j)%subshell_two_s        = term_ls_i2(i)%subshell_two_s
                end do
            case default
                stop  "enumerate_ls_subshell_terms(): program stop C."
            end select
        case (7)
            select case (n)
            case (1)
                i = 1;       j = 1
                ls(j)%subshell_l  = term_ls_k1(i)%subshell_l
                ls(j)%w        = term_ls_k1(i)%w
                ls(j)%subshell_two_quasispin        = term_ls_k1(i)%subshell_two_quasispin
                ls(j)%subshell_two_l       = term_ls_k1(i)%subshell_two_l
                ls(j)%subshell_two_s        = term_ls_k1(i)%subshell_two_s
            case (2)
                do i = 1,15
                    j = j + 1
                    ls(j)%subshell_l  = term_ls_k2(i)%subshell_l
                    ls(j)%w        = term_ls_k2(i)%w
                    ls(j)%subshell_two_quasispin        = term_ls_k2(i)%subshell_two_quasispin
                    ls(j)%subshell_two_l       = term_ls_k2(i)%subshell_two_l
                    ls(j)%subshell_two_s        = term_ls_k2(i)%subshell_two_s
                end do
            case default
                stop  "enumerate_ls_subshell_terms(): program stop D."
            end select
        case (8)
            select case (n)
            case (1)
                i = 1;       j = 1
                ls(j)%subshell_l  = term_ls_l1(i)%subshell_l
                ls(j)%w        = term_ls_l1(i)%w
                ls(j)%subshell_two_quasispin        = term_ls_l1(i)%subshell_two_quasispin
                ls(j)%subshell_two_l       = term_ls_l1(i)%subshell_two_l
                ls(j)%subshell_two_s        = term_ls_l1(i)%subshell_two_s
            case (2)
                do i = 1,17
                    j = j + 1
                    ls(j)%subshell_l  = term_ls_l2(i)%subshell_l
                    ls(j)%w        = term_ls_l2(i)%w
                    ls(j)%subshell_two_quasispin        = term_ls_l2(i)%subshell_two_quasispin
                    ls(j)%subshell_two_l       = term_ls_l2(i)%subshell_two_l
                    ls(j)%subshell_two_s        = term_ls_l2(i)%subshell_two_s
                end do
            case default
                stop  "enumerate_ls_subshell_terms(): program stop E."
            end select
        case (9,10)
            ! Equivalent electrons give singlets at even L and triplets at odd L.
            select case (n)
            case (1)
                j = 1
                ls(j) = subshell_term_ls(l_shell,0,2*l_shell,2*l_shell,1)
            case (2)
                j = 2*l_shell+1
                ls(1) = subshell_term_ls(l_shell,0,2*l_shell+1,0,0)
                do i = 1, 2*l_shell
                    ls(i+1) = subshell_term_ls(l_shell,0,2*l_shell-1,2*i,2*mod(i,2))
                end do
            case default
                error stop "JJ2LSJ m and n shells require one or two electrons"
            end select
        case default
            stop  "enumerate_ls_subshell_terms(): program stop G."
        end select
        number = j
    end subroutine enumerate_ls_subshell_terms

    ! The input from the screen.
    !
    ! inscreenlev.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: Dec 2015
    ! Modified by G. Gaigalas. 2022
    ! Modified by G. Gaigalas. 2024
    subroutine read_ls_conversion_options(thresh,levels,number_of_levels,write_csf_list,write_j_vectors,  &
        unique,transformation_output_mode)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real32
        use orbital_csf_basis,           only: csf_count, orbital_count
        use selected_level_indices,          only: active_level_count
        use standard_io_units, only: stderr_unit
        use numeric_constants,          only: angular_tolerance, zero
        use transition_basis_blocks,           only: transition_level_counts_by_block, transition_csf_counts_by_block, &
            transition_block_count, transition_two_j_by_block
        use csf_pair_angular_state,             only: core_orbital_count
        use calculation_state_files, only: read_state_name, lsj_label_path
        use atomic_calculation_state,           only: nuclear_charge, electron_count

        ! Procedure interfaces
        use read_mixing_blocks_interface, only: read_mixing_blocks

        ! Module dependencies
        use yes_no_input, only: read_yes_no
        use file_open, only: open_calculation_file
        use csf_full_basis_reader, only: load_csf_basis
        implicit none (type, external)

        ! Arguments and local declarations
        save  ! Preserve the legacy serial controller workspace.
        integer, intent(out)      :: write_csf_list,write_j_vectors,unique,transformation_output_mode
        real(real64), intent(out) :: thresh
        integer, intent(out) :: number_of_levels(max_symmetry_blocks)
        integer, intent(out) :: levels(max_symmetry_blocks, max_selected_levels)
        real(real64) :: phase_start
        integer            :: i, ii, isum, k, nci, io_status
        integer            :: iblock, number_of_levels_tmp
        logical            :: yes, fail
        character(len=24)  :: name
        character(len=256) :: record, util_csl_file
        integer    :: posi(max_symmetry_blocks)
        integer :: levels_tmp(1:max_selected_levels)
        transition_block_count = 0
        call read_state_name(name)
        k = len_trim(name) + 1
        ! Load the configuration state functions from the descriptor pair.
        phase_start=profile_time()
        call load_csf_basis(name,core_orbital_count,io_status)
        if (io_status /= 0) error stop "JJ2LSJ descriptor basis load failed"
        call profile_emit("csf_read_parse",phase_start,items=csf_count)
        phase_start=profile_time()
        allocate(intermediate_two_j(1:orbital_count,1:csf_count))
        call select_ls_transition_labels
        call profile_emit("jcoup",phase_start,items=csf_count)
        write (stderr_unit,*)
        write (stderr_unit,*) "Mixing coefficients from a CI calc.?"
        yes = read_yes_no ()
        if (yes) then
            nci = 0
        else
            nci = 1
        end if
        ! 2017
        write (stderr_unit,*) "Do you need a unique labeling? (y/n)"
        yes = read_yes_no ()
        if (yes) then
            unique = 1
        else
            unique = 0
        end if
        ! 2017 end
        !
        !    Get the eigenvectors
        phase_start=profile_time()
        call read_mixing_blocks (name, nci)
        call profile_emit("mix_read",phase_start,items=active_level_count)
        write (stderr_unit,*) "Default settings?  (y/n) "
        yes = read_yes_no ()
        if (yes) then
            coefficient_threshold = 0.005_real64
            isum = 0
            limit_omitted_weight = 1
            max_omitted_percent = 1
            thresh = 0.001_real64
            write_csf_list = 0
            write_j_vectors = 0
            transformation_output_mode = 0
            do i = 1, transition_block_count
                number_of_levels(i) = transition_level_counts_by_block(i)
                if(transition_level_counts_by_block(i) /= 0) then
                    do ii = 1, transition_level_counts_by_block(i)
                        isum = isum + 1
                        levels(i,ii) = isum
                    end do
                end if
            end do
            ! 2024
            ls_convention = 1
        ! 2024 end
        else
            write (stderr_unit,*) "All levels (Y/N)"
            yes = read_yes_no ()
            write (stderr_unit,*)
            if (yes) then
                isum = 0
                do i = 1, transition_block_count
                    number_of_levels(i) = transition_level_counts_by_block(i)
                    if(transition_level_counts_by_block(i) /= 0) then
                        do ii = 1, transition_level_counts_by_block(i)
                            isum = isum + 1
                            levels(i,ii) = isum
                        end do
                    end if
                end do
            else
                number_of_levels = 0
                posi(1) = 0
                do i = 1, transition_block_count-1
                    posi(i+1) = transition_level_counts_by_block(i) + posi(i)
                end do
                write (stderr_unit,*) "Maximum number of ASFs is:",active_level_count
                write (stderr_unit,*) "Enter the level numbers of the ASF which are to be transformed,"
                do
                    write (stderr_unit,*) "Enter the block number"
                    read (*, "(I2)") iblock
                    write (stderr_unit,*) "The block number is:",iblock
                    write (stderr_unit,*) " e.g. 1 3 4  7 - 20  48  69 - 85 :"
                    read (*, "(a)") record
                    call parse_level_selection(record,levels_tmp,number_of_levels_tmp,fail)
                    if (fail) then
                        write (stderr_unit,*) "Unable to interprete the serial level numbers; redo ..."
                        cycle
                    end if
                    number_of_levels(iblock) = number_of_levels_tmp
                    if (active_level_count < number_of_levels(iblock)) then
                        write (stderr_unit,*) "There are to much ASF:", number_of_levels(iblock)
                        cycle
                    end if
                    do i = 1,number_of_levels_tmp
                        levels(iblock,i) = posi(iblock) + levels_tmp(i)
                    end do
                    write (stderr_unit,*) ""
                    write (stderr_unit,*) "Do you need to include more levels?  (y/n)"
                    yes = read_yes_no ()
                    if (.not. yes) exit
                end do
            end if
            do
                write (stderr_unit,*) "Maximum % of omitted composition"
                ! 2024
                write (stderr_unit,*) "(default setting: 1.000)"
                ! 2024 end
                read *, max_omitted_percent
                if(max_omitted_percent == zero) then
                    limit_omitted_weight = 0
                    coefficient_threshold = angular_tolerance*angular_tolerance
                    write_csf_list = 1
                    write_j_vectors = 1
                else if( max_omitted_percent > zero) then
                    limit_omitted_weight = 1
                    write (stderr_unit,*) "What is the value below which an eigenvector component"
                    write (stderr_unit,*) "is to be neglected in the determination of the LSJ expansion:"
                    write (stderr_unit,"(A,F8.5)") " should be smaller than:",max_omitted_percent*0.01_real64
                    ! 2024
                    write (stderr_unit,*) "(default setting: 0.005)"
                    ! 2024 end
                    read *, coefficient_threshold
                else
                    write (stderr_unit,*) " The maximum of omitted composition can be 100%"
                    cycle
                end if
                exit
            end do
            write (stderr_unit,*) "What is the value below which an eigenvector composition"
            write (stderr_unit,*) "is to be neglected for printing?"
            ! 2024
            write (stderr_unit,*) "(default setting: 0.001)"
            ! 2024 end
            read *, thresh
            if( max_omitted_percent > zero) then
                write (stderr_unit,*) "Do you need the output file *.lsj.c?  (y/n)"
                yes = read_yes_no ()
                if (yes) then
                    write_csf_list = 1
                else
                    write_csf_list = 0
                end if
                write (stderr_unit,*) "Do you need the output file *.lsj.j?  (y/n)"
                yes = read_yes_no ()
                if (yes) then
                    write_j_vectors = 1
                else
                    write_j_vectors = 0
                end if
            end if
            if(max_omitted_percent == zero) then
                write (stderr_unit,*)                                            &
                    "Do you need the transformation output file *.lsj.transform.cache?  (y/n)"
                yes = read_yes_no ()
                if (yes) then
                    transformation_output_mode = 1
                else
                    write (stderr_unit,*)                                         &
                        "Will you use the transformation file *.lsj.transform.cache?  (y/n)"
                    yes = read_yes_no ()
                    if (yes) then
                        transformation_output_mode = 2
                    else
                        transformation_output_mode = 0
                    end if
                end if
            else
                transformation_output_mode = 0
            end if
            ! 2024
            write (stderr_unit,*)                                               &
                "Which convention should be used for d shell identification?"
            write (stderr_unit,*)                                               &
                "Seniority notation / Nielson-Koster indexes (0/1)"
            read *, ls_convention
            ! 2024 end
        end if
        write (stderr_unit,*)
        if (limit_omitted_weight == 1)             &
            write (*,"(A,F8.3)") " Maximum % of omitted composition is ",max_omitted_percent
        write (*,"(A,ES8.1,A)")   &
            " Below ",coefficient_threshold," the eigenvector component is to be neglected for calculating"
        write (*,"(A,ES8.1,A)")   &
            " Below ",thresh," the eigenvector composition is to be neglected for printing"
        ! 2024
        write (stderr_unit,*)
        if (ls_convention == 0) then
            write (*,"(A)")                                                &
                " Seniority notaton is used for d shell identification in"
            write (*,"(5A)")"   ",lsj_label_path(name, nci == 0), ", ",       &
                name(1:k-1)//".uni.lsj"//".lbl",","
            write (*,"(5A)")"   ",name(1:k-1)//".uni.lsj"//".sum",", and ",&
                name(1:k-1)//".lsj"//".j","  files"
        else
            write (*,"(A)")                                                &
                " Nielson-Koster indexes is used for d shell identification in"
            write (*,"(5A)")"   ",lsj_label_path(name, nci == 0), ", ",       &
                name(1:k-1)//".uni.lsj"//".lbl",","
            write (*,"(5A)")"   ",name(1:k-1)//".uni.lsj"//".sum",", and ",&
                name(1:k-1)//".lsj"//".j","  files"
        end if
        write (*,"(A)")                                                  &
            " while seniority notation for this shell is used in other files"
        ! 2024 end
        print *, " "
        if (transformation_output_mode /= 0) call open_transformation_cache(name, transformation_output_mode)
        ! Opening the files  *.lsj.c
        if(write_csf_list == 1) then
            util_csl_file = name(1:k-1)//".lsj"//".c"
            call open_calculation_file(ls_csf_output_unit,util_csl_file,"formatted","UNKNOWN",io_status)
            if (io_status == 1) then
                print *, "Error when opening",util_csl_file
                stop
            end if
        end if
        ! Open the source-specific LSJ label report.
        util_csl_file = lsj_label_path(name, nci == 0)
        call open_calculation_file(lsj_label_unit,util_csl_file,"formatted","UNKNOWN",io_status)
        if (io_status == 1) then
            print *, "Error when opening",util_csl_file
            stop
        end if
        ! 2017
        if (unique == 1) then
            util_csl_file = name(1:k-1)//".uni.lsj"//".lbl"
            call open_calculation_file(unique_label_unit,util_csl_file,"formatted","UNKNOWN",io_status)
            if (io_status == 1) then
                print *, "Error when opening",util_csl_file
                stop
            end if
            util_csl_file = name(1:k-1)//".uni.lsj"//".sum"
            call open_calculation_file(unique_summary_unit,util_csl_file,"formatted","UNKNOWN",io_status)
            if (io_status == 1) then
                print *, "Error when opening",util_csl_file
                stop
            end if
        end if
        ! 2017 end
        write(lsj_label_unit,"(A53)") " Pos   J   Parity      Energy Total      Comp. of ASF"
        ! Opening the files   *.lsj.j and
        if(write_j_vectors == 1) then
            util_csl_file = name(1:k-1)//".lsj"//".j"
            call open_calculation_file(j_vector_output_unit,util_csl_file,"formatted","UNKNOWN",io_status)
            if (io_status == 1) then
                print *, "Error when opening",util_csl_file
                stop
            end if
            write (j_vector_output_unit,"(2X,A6,A,F5.1,A,I3,A,I7)" )    &
                name(1:k-1),"  Z = ",nuclear_charge ,"  NEL = ",electron_count,"   NCFG ="! ,asf_set_LS%csf_set_LS%nocsf
        end if
    end subroutine read_ls_conversion_options

    subroutine open_transformation_cache(name, mode)
        use atomic_calculation_state, only: electron_count
        use csf_descriptor_reader, only: loaded_csf_basis_id
        use orbital_csf_basis, only: csf_count, orbital_count
        use transition_basis_blocks, only: transition_block_count, transition_csf_counts_by_block, &
            transition_level_counts_by_block, transition_two_j_by_block

        character(*), intent(in) :: name
        integer, intent(in) :: mode
        character(:), allocatable :: path
        character(len(transformation_signature)) :: signature
        character(len(loaded_csf_basis_id)) :: basis_id
        integer :: status, block, convention, dimensions(4), block_header(4), stored(4)

        path = trim(name)//".lsj.transform.cache"
        dimensions = [electron_count, csf_count, orbital_count, transition_block_count]
        if (mode == 1) then
            open(newunit=transformation_unit, file=path, form="unformatted", status="new", action="write", iostat=status)
            if (status /= 0) error stop "Cannot create JJ2LSJ transformation cache: "//path
            write(transformation_unit) transformation_signature, loaded_csf_basis_id, ls_convention, dimensions
        else
            open(newunit=transformation_unit, file=path, form="unformatted", status="old", action="read", iostat=status)
            if (status /= 0) error stop "Cannot open JJ2LSJ transformation cache: "//path
            read(transformation_unit, iostat=status) signature, basis_id, convention, stored
            if (status /= 0) error stop "Invalid JJ2LSJ transformation cache header: "//path
            if (signature /= transformation_signature) error stop "Unsupported JJ2LSJ transformation cache: "//path
            if (basis_id /= loaded_csf_basis_id) error stop "JJ2LSJ transformation cache CSF basis differs: "//path
            if (convention /= ls_convention .or. any(stored /= dimensions)) &
                error stop "JJ2LSJ transformation cache convention or dimensions differ: "//path
        end if
        do block = 1, transition_block_count
            block_header = [block, transition_csf_counts_by_block(block), transition_level_counts_by_block(block), &
                transition_two_j_by_block(block)+1]
            if (mode == 1) then
                write(transformation_unit) block_header
            else
                read(transformation_unit, iostat=status) stored
                if (status /= 0) error stop "Incomplete JJ2LSJ transformation cache directory: "//path
                if (any(stored /= block_header)) error stop "JJ2LSJ transformation cache block metadata differ: "//path
            end if
        end do
    end subroutine open_transformation_cache

    subroutine close_transformation_cache(mode)
        integer, intent(in) :: mode
        integer :: status
        character(len(transformation_done)) :: completion

        if (mode == 1) then
            write(transformation_unit) transformation_done
        else
            read(transformation_unit, iostat=status) completion
            if (status /= 0) error stop "Missing JJ2LSJ transformation cache completion record"
            if (completion /= transformation_done) error stop "Unexpected records in JJ2LSJ transformation cache"
        end if
        close(transformation_unit)
    end subroutine close_transformation_cache

    !     Attempts to interprete the serial level numbers which are given
    !     in record. These level numbers can be given in the format:
    !                1 3 4  7 - 20  48  69 - 85
    !     Any order and 'overlapping' intervals are also supported.
    !     The procedure returns with fail = .true. if the level numbers
    !     cannot be interpreted properly (fail = .false. otherwise).
    !
    !     The level numbers are returned in the vector
    !                                          levels(1:number_of_levels). *
    !     The procedure assumes that this vector has a sufficient
    !     dimension to store all level numbers.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: Dec 2015
    subroutine parse_level_selection(record,levels,number_of_levels,fail)
        ! Procedure interfaces
        use ascii_digit_value_interface, only: ascii_digit_value
        implicit none (type, external)

        ! Arguments and local declarations
        save  ! Preserve the legacy serial controller workspace.
        character(len=*), intent(in) :: record
        logical, intent(out)         :: fail
        integer, intent(out)         :: number_of_levels
        integer, intent(inout) :: levels(1:max_selected_levels)
        logical   :: low_to(2000)
        character(len=500)        :: string
        integer                   :: a, i, lower, n
        integer   :: low(2000)
        integer :: run(max_selected_levels)
        fail = .true.;  levels(1:max_selected_levels) = 0;  number_of_levels = 0
        run(1:max_selected_levels) = 0;   string = adjustl(record)
        n = 0;   lower = 0;   low_to(:) = .false.
        do
            string = adjustl(string)
            if (string(1:1) == " ") then
                if (n == 0) return
                if (low_to(n)) return
                exit
            else if (string(1:1) == "-") then
                if (n == 0)  return
                low_to(n)   = .true.;    string(1:1) = " "
            else if (string(1:1) == "0"   .or.   string(1:1) == "1" .or.  &
                string(1:1) == "2"   .or.   string(1:1) == "3" .or.  &
                string(1:1) == "4"   .or.   string(1:1) == "5" .or.  &
                string(1:1) == "6"   .or.   string(1:1) == "7" .or.  &
                string(1:1) == "8"   .or.   string(1:1) == "9") then
                a = ascii_digit_value(string(1:1));    lower = 10*lower + a
                if (string(2:2) == " "   .or.   string(2:2) == "-") then
                    n = n + 1;    low(n) = lower;    lower = 0
                end if
                string(1:1) = " "
            else
                exit
            end if
        end do
        ! Determine no_eigenpairs and max_eigenpair
        levels(:) = 0
        do  i = 1,n
            if (low_to(i)) then
                if (low(i) <= low(i+1)) then;   run(low(i):low(i+1)) = 1
                else;                           run(low(i+1):low(i)) = 1
                end if
                if (low_to(i+1)) then;   return
                else;                    cycle
                end if
            else
                run(low(i)) = 1
            end if
        end do
        number_of_levels = 0
        do  i = 1,max_selected_levels
            if (run(i) == 1) then
                number_of_levels         = number_of_levels + 1
                levels(number_of_levels) = i
            end if
        end do
        fail = .false.
    end subroutine parse_level_selection

    ! Controls the transformation of atomic states from a jj-
    ! to a LS-coupled CSF basis.
    !
    !        setLS.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: Dec 2015
    ! Modified by G. Gaigalas,. May 2021
    ! Modified by C. Cychen. 2021
    ! Modified by G. Gaigalas. 2022
    subroutine run_ls_coupling_conversion
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real32
        use, intrinsic :: iso_fortran_env, only: int64
        use eigenstate_expansions,          only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
        use orbital_csf_basis,           only: csf_count
        use selected_level_indices,          only: active_level_indices, active_level_count
        use numeric_constants, only: zero
        use transition_basis_blocks, only: transition_csf_counts_by_block, transition_block_count
        use csf_quantum_queries, only: csf_total_j_code, csf_parity_code
        use half_integer_format, only: format_half_integer
        implicit none (type, external)

        ! Local declarations
        save  ! Preserve the legacy serial controller workspace.
        ! GG NIST
        integer :: i, j, jj, string_l,symmetry_block_index,write_csf_list,write_j_vectors,unique,transformation_output_mode
        integer :: level, nocsf_min, lev, string_length
        integer :: nocsf_max, sum_nocsf_min, before_j
        ! GG NIST
        integer :: first_csf_index, last_csf_index, ncf_ls_jj_max
        type(ls_ranking), allocatable :: rankings(:)
        integer :: post_lev, post_workers
        real(real64) :: rank_start
        real(real64) :: phase_start, format_start, format_seconds
        real(real64) :: thresh, wa, wb,ssms,g_j,g_jls, sum_absolute_coefficients
        character(len=4)   :: string_cnum
        character(len=164)  :: string_csf_one
        integer      :: number_of_levels(max_symmetry_blocks)
        integer :: levels(max_symmetry_blocks, max_selected_levels)
        integer     :: max_comp(1:max_selected_levels)
        integer, pointer     :: leading_ls(:)
        integer, pointer       :: iw(:)
        real(real64), pointer  :: weights(:), weights2(:)
        integer, pointer  :: retained_csf_mask(:)
        character(len=164) :: string_csf(1:max_selected_levels)
        integer(int64) :: level_mixing_offset, csf_count_int64
        ssms = zero;   g_j = zero;   g_jls = zero;    before_j = 0
        call read_ls_conversion_options(thresh,levels,number_of_levels,write_csf_list,write_j_vectors,unique, &
            transformation_output_mode)
        allocate(retained_csf_mask(csf_count))
        csf_count_int64 = csf_count
        do  symmetry_block_index = 1, transition_block_count
            if(symmetry_block_index == 1) then
                first_csf_index = 1
                last_csf_index = transition_csf_counts_by_block(symmetry_block_index)
            else
                first_csf_index = last_csf_index + 1
                last_csf_index = first_csf_index + transition_csf_counts_by_block(symmetry_block_index) - 1
            end if
            if (number_of_levels(symmetry_block_index) == 0) cycle
            phase_start=profile_time()
            retained_csf_mask = 0
            do  i = first_csf_index, last_csf_index
                sum_absolute_coefficients = zero
                do  lev = 1, number_of_levels(symmetry_block_index)
                    level = levels(symmetry_block_index,lev)
                    level_mixing_offset = (level-1)*csf_count_int64
                    sum_absolute_coefficients = sum_absolute_coefficients + dabs(csf_mixing_coefficients(i+level_mixing_offset))
                end do
                if(dabs(sum_absolute_coefficients) >= dabs(coefficient_threshold)) retained_csf_mask(i) = 1
            end do
            call profile_emit("screen",phase_start,symmetry_block_index,items= &
                & count(retained_csf_mask(first_csf_index:last_csf_index)==1))
            phase_start=profile_time()
            write(*,*)"The program generates a list of CSFs in LS-coupling"
            call build_ls_basis(retained_csf_mask,first_csf_index,last_csf_index)
            call profile_emit("ls_basis",phase_start,symmetry_block_index,items=asf_set_ls%csf_set_ls%ls_csf_count)
            phase_start=profile_time()
            ! output to *.lsj.c
            if(write_csf_list == 1) call write_ls_csf_list (ls_csf_output_unit,asf_set_ls%csf_set_ls,before_j)
            call profile_emit("csf_write",phase_start,symmetry_block_index)
            phase_start=profile_time()
            allocate(asf_set_ls%asf(1:active_level_count))
            do i = 1, active_level_count
                allocate(asf_set_ls%asf(i)%eigenvector(1:asf_set_ls%csf_set_ls%ls_csf_count))
            end do
            do lev = 1, number_of_levels(symmetry_block_index)
                level = levels(symmetry_block_index,lev)
                asf_set_ls%asf(level)%level_no = level
                asf_set_ls%asf(level)%energy = hamiltonian_energy_shift+shifted_eigenvalues(level)
                asf_set_ls%asf(level)%eigenvector = zero
            end do
            asf_set_ls%asf_count = number_of_levels(symmetry_block_index)
            nocsf_max = asf_set_ls%csf_set_ls%ls_csf_count
            ! GG NIST
            ncf_ls_jj_max = max(1000,last_csf_index-first_csf_index+1,asf_set_ls%csf_set_ls%ls_csf_count+1)
            allocate (weights(ncf_ls_jj_max))
            allocate (weights2(ncf_ls_jj_max))
            allocate (iw(ncf_ls_jj_max))
            allocate (leading_ls(ncf_ls_jj_max))
            csf_count_int64 = csf_count
            call profile_emit("workspace",phase_start,symmetry_block_index)
            do  lev = 1, number_of_levels(symmetry_block_index)
                level = levels(symmetry_block_index,lev)
                phase_start=profile_time()
                weights = zero;    iw = 0;    wb = zero
                do  i = first_csf_index, last_csf_index
                    if(retained_csf_mask(i) == 1) then
                        level_mixing_offset = (level-1)*csf_count_int64
                        wa = csf_mixing_coefficients(i+level_mixing_offset) * csf_mixing_coefficients(i+level_mixing_offset)
                        wb = wb + wa
                        do  j = 1,999
                            if (wa > weights(j)) then
                                weights(j+1:1000) = weights(j:999)
                                weights(j)       = wa
                                iw(j+1:1000)      = iw(j:999)
                                iw(j)            = i - first_csf_index + 1
                                exit
                            end if
                        end do
                    end if
                end do
                call profile_emit("jj_rank",phase_start,symmetry_block_index,level)
                if(lev == 1) then
                    print *, " "
                    write(*,"(A,A)") " .  .  .  .  .  .  .  .  .  .  .  .  .  .  .",&
                        "  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  ."
                    write(*, "(A,2X,I4,16X,A,I3)")                            &
                        " Under investigation is the block:",symmetry_block_index,             &
                        " The number of eigenvectors:", number_of_levels(symmetry_block_index)
                    write(*,"(A,I10,10X,A,I10)")                              &
                        " The number of CSF (in jj-coupling):",transition_csf_counts_by_block(symmetry_block_index), &
                        " The number of CSF (in LS-coupling):",asf_set_ls%csf_set_ls%ls_csf_count
                else
                    print *, " "
                    write(*,"(A)") " .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  ."
                    write(*,"(A)") " The new level is under investigation."
                end if
                ! perform the transformation
                if(lev == 1) then
                    write(*,*)                                              &
                        "The program generates a list of transformation matrix"
                    phase_start=profile_time()
                    call transform_asf_to_ls                                             &
                        (iw(1),retained_csf_mask,number_of_levels(symmetry_block_index),symmetry_block_index,levels,&
                        first_csf_index,last_csf_index,transformation_output_mode)
                    call profile_emit("transform",phase_start,symmetry_block_index)
                    allocate(rankings(number_of_levels(symmetry_block_index)))
                    post_workers=post_worker_count(number_of_levels(symmetry_block_index))
                    phase_start=profile_time()
!$omp parallel do default(none) schedule(static) num_threads(post_workers) if(post_workers>1) &
!$omp shared(rankings,number_of_levels,symmetry_block_index,levels,asf_set_LS) private(post_lev,rank_start)
                    do post_lev=1,number_of_levels(symmetry_block_index)
                        rank_start=profile_time()
                        call rank_ls_components(asf_set_ls%asf(levels(symmetry_block_index,post_lev))%eigenvector, &
                            & rankings(post_lev))
                        rankings(post_lev)%seconds=profile_time()-rank_start
                    end do
!$omp end parallel do
                    call profile_emit("ls_rank_wall",phase_start,symmetry_block_index,items=post_workers)
                    do post_lev=1,number_of_levels(symmetry_block_index)
                        call profile_emit("ls_rank_asf",phase_start,symmetry_block_index,levels(symmetry_block_index,post_lev), &
                            asf_set_ls%csf_set_ls%ls_csf_count,elapsed=rankings(post_lev)%seconds)
                    end do
                end if
                ! output to the screen jj- coupling
                phase_start=profile_time()
                print *, "Weights of major contributors to ASF in jj-coupling:"
                print *, " "
                print *, " Level  J Parity      CSF contributions"
                print *, " "
                if (level > 1   .and.   dabs(wb) > 1.0001_real64) then
                    print *, "level, wb = ",level,wb
                    stop "run_ls_coupling_conversion(): program stop A."
                end if
                nocsf_min = 5
                do j = 1,5
                    if(abs(weights(j)) < 0.00001_real64) then
                        nocsf_min = j - 1
                        exit
                    end if
                end do
                ! CYC-2023 Corrections, fix the wrong output information.
                if(csf_parity_code(iw(1)+first_csf_index-1) == -1) then
                    asf_set_ls%asf(level)%parity = "-"
                else if(csf_parity_code(iw(1)+first_csf_index-1) == 1) then
                    ! CYC-2023 end
                    asf_set_ls%asf(level)%parity = "+"
                else
                    stop "run_ls_coupling_conversion: program stop D."
                end if
                asf_set_ls%asf(level)%two_j = csf_total_j_code(iw(1)+first_csf_index-1) - 1
                call format_half_integer(asf_set_ls%asf(level)%two_j,string_cnum,string_l)
                print '(1x,i4,1x,2a4,2x,100(3x,f8.5," of",i7,3x,f8.5," of",i7,3x,f8.5," of",i7,3x,f8.5," of",i7/16x))', &
                    active_level_indices(level),string_cnum(1:string_l),               &
                    asf_set_ls%asf(level)%parity,(weights(j),iw(j),j=1,nocsf_min)
                print*, "              Total sum over  weight (in jj) is:",wb
                print *, " "
                print *, "Definition of leading CSF:"
                print *, " "
                ! CGG            call prCSFjj(-1,iw(1))
                call write_jj_csf(-1,iw(1)+first_csf_index-1,iw(1))
                ! output to the screen LS- coupling
                print *, " "
                print *, " "
                print *, "Weights of major contributors to ASF in LS-coupling:"
                print *, " "
                print *, " Level  J Parity      CSF contributions"
                print *, " "
                call profile_emit("jj_stdout",phase_start,symmetry_block_index,level)
                phase_start=profile_time()
                sum_nocsf_min = 0
                weights = zero;  weights2 = zero;  iw = 0;  wb = zero
                wb=rankings(lev)%norm
                do i=1,asf_set_ls%csf_set_ls%ls_csf_count
                    iw(i)=rankings(lev)%order(i)
                    weights(i)=asf_set_ls%asf(level)%eigenvector(iw(i))
                    weights2(i)=weights(i)*weights(i)
                end do
                call profile_emit("ls_rank_gather",phase_start,symmetry_block_index,level,asf_set_ls%csf_set_ls%ls_csf_count)
                phase_start=profile_time()
                if(limit_omitted_weight == 1) then
                    if (level > 1 .and. (wb*100.0_real32+dabs(max_omitted_percent)) <= 100.0_real32) then
                        write (*,"(A,F5.0,A)")   &
                            "Attention!!!  The program not reach the accuracy",max_omitted_percent," %"
                        print *, "level, sum of weights in LSJ = ",level,wb
                        print *,""
                    end if
                else
                    if (level > 1 .and. abs(wb) > 1.001_real64) then
                        print*,"Attention!!!  The sum of weights is bigger then 1.0001"
                        print *, "level, sum of weights in LSJ = ",level,wb
                    end if
                    if (level > 1 .and.  abs(wb) < 0.998_real64) then
                        print*,"Attention!!!  The sum of weights is little then .998"
                        print *, "level, sum of weights in LSJ = ",level,wb
                    end if
                end if
                nocsf_min = nocsf_max;   jj = 0
                do j = 1,nocsf_max
                    if(dabs(weights2(j)) < dabs(thresh)) then
                        nocsf_min = j - 1
                        exit
                    end if
                    jj = jj + 1
                    leading_ls(jj) = iw(j)
                end do
                max_comp(level) = iw(1)
                sum_nocsf_min = sum_nocsf_min + nocsf_min
                if (nocsf_min <= 0) then
                    stop "run_ls_coupling_conversion(): program stop C."
                end if
                call profile_emit("candidate_select",phase_start,symmetry_block_index,level,nocsf_min)
                phase_start=profile_time()
                call format_half_integer(asf_set_ls%csf_set_ls%csf(iw(1))%two_j,string_cnum,string_l)
                print '(1x,i4,1x,2a4,2x,100(3x,f8.5," of",i7,3x,f8.5," of",i7,3x,f8.5," of",i7,3x,f8.5," of",i7/16x))', &
                    asf_set_ls%asf(level)%level_no,string_cnum(1:string_l), &
                    asf_set_ls%csf_set_ls%csf(iw(1))%parity,                   &
                    (weights2(j),iw(j),j=1,5)
                print*, "              Total sum over  weight (in LSJ) is:",wb
                print *, " "
                print *, "Definition of leading CSF:"
                print *, " "
                if(sum_nocsf_min > 5) sum_nocsf_min = 5
                do  i = 1,csf_count
                    do  j = 1, sum_nocsf_min
                        if (i == iw(j)) then
                            call write_ls_csf (-1,asf_set_ls%csf_set_ls,leading_ls(j))
                            exit
                        end if
                    end do
                end do
                ! Output the LSJ label report.
                call profile_emit("ls_stdout",phase_start,symmetry_block_index,level)
                phase_start=profile_time()
                format_seconds=0.0_real64
                if(before_j == 0) then
                    before_j = 1
                else if(lev == 1 .and. before_j == 1) then
                    write(lsj_label_unit,"(A1)") " "
                end if
                do j = 1,nocsf_min
                    format_start=profile_time()
                    call pack_ls_csf_label(asf_set_ls%csf_set_ls,iw(j),string_csf_one)
                    format_seconds=format_seconds+(profile_time()-format_start)
                    if(j == 1) then
                        write(lsj_label_unit,"(I3,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)")        &
                            asf_set_ls%asf(level)%level_no,string_cnum(1:string_l),&
                            asf_set_ls%csf_set_ls%csf(iw(1))%parity,               &
                            asf_set_ls%asf(level)%energy,wb*100,"%"
                        string_csf(level) = string_csf_one
                    end if
                    string_length = len_trim(string_csf_one)
                    write(lsj_label_unit,"(7X,F12.8,3X,F11.8,3X,A)") weights(j),weights2(j),string_csf_one(1:string_length)
                end do
                ! output  *.lsj.j and
                call profile_emit("label_format",phase_start,symmetry_block_index,level,nocsf_min,elapsed=format_seconds)
                call profile_emit("label_write",phase_start,symmetry_block_index,level, &
                    elapsed=profile_time()-phase_start-format_seconds)
                phase_start=profile_time()
                if(write_j_vectors == 1) then
                    if(lev == 1) then
                        write (j_vector_output_unit, "(//A8,I4,2X,A8,I4)" ) "  2*J = ",       &
                            asf_set_ls%asf(level)%two_j,"NUMBER =",number_of_levels(symmetry_block_index)
                    end if
                    string_length = len_trim(string_csf(level))
                    write(j_vector_output_unit,"(3(A8,F15.10))") "Ssms=", ssms, "g_J=",g_j, "g_JLS=", g_jls
                    write(j_vector_output_unit,"(I6,F16.9,2X,A)")                                &
                        max_comp(level),asf_set_ls%asf(level)%energy,string_csf(level)(1:string_length)
                    write(j_vector_output_unit,"(7F11.8)")                                       &
                        (asf_set_ls%asf(level)%eigenvector(i), i=1,asf_set_ls%csf_set_ls%ls_csf_count)
                end if
                call profile_emit("vector_write",phase_start,symmetry_block_index,level)
            end do
            phase_start=profile_time()
            call release_ls_expansion(asf_set_ls)
            deallocate(weights)
            deallocate(weights2)
            deallocate(iw)
            deallocate(leading_ls)
            deallocate(rankings)
            call profile_emit("release",phase_start,symmetry_block_index)
        end do
        if (transformation_output_mode /= 0) call close_transformation_cache(transformation_output_mode)
        deallocate(retained_csf_mask)
        if(write_csf_list == 1) call profile_file_size("csf_bytes",ls_csf_output_unit)
        if (ls_csf_output_unit /= -1) close(ls_csf_output_unit)
        ls_csf_output_unit = -1
        ! 2017
        if (unique == 1) then
            phase_start=profile_time()
            call assign_unique_lsj_labels(lsj_label_unit, unique_summary_unit, unique_label_unit)
            call profile_emit("unique_total",phase_start)
            call profile_file_size("unique_summary_bytes",unique_summary_unit)
            call profile_file_size("unique_label_bytes",unique_label_unit)
            close(unique_summary_unit)
            close(unique_label_unit)
        end if
        ! 2017
        call profile_file_size("label_bytes",lsj_label_unit)
        close(lsj_label_unit)
        if(write_j_vectors == 1) write(j_vector_output_unit,"(A)") "**"
        if(write_j_vectors == 1) call profile_file_size("vector_bytes",j_vector_output_unit)
        if (j_vector_output_unit /= -1) close(j_vector_output_unit)
        j_vector_output_unit = -1
        deallocate(intermediate_two_j)
        return
    end subroutine run_ls_coupling_conversion

    !     Encording all CSF in LS - coupling.
    !
    !     Rules for encoding
    !        1. All blanks deleted
    !        2. If Qi=1, omit Qi
    !        3. If Qi=1 or Qi>=4l+1, omit ALFAi
    !        4. If i=1 or (Qi=4l+2 and i<>m), insert '.'; else _BETAi.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: Dec 2015
    subroutine pack_ls_csf_label(csf_set_ls,csf_index,string_csf)
        ! Procedure interfaces
        use pack_ls_subshell_label_interface, only: pack_ls_subshell_label
        implicit none (type, external)

        ! Arguments and local declarations
        save  ! Preserve the legacy serial controller workspace.
        type(csf_basis_ls), intent(in) :: csf_set_ls
        integer, intent(in)            :: csf_index
        character(len=164), intent(out) :: string_csf
        integer                                 :: i, j, counter
        integer, pointer          :: occupation(:)
        integer                   :: q(18)
        character(len=4)                        :: ls, xls
        character(len=4),pointer  :: string_ls(:), string_xls(:)
        character(len=3)         :: couple(35)
        character(len=3)          :: elc(18)
        character(len=2)                        :: string
        character(len=164)                       :: string_csf_pack
        character(len=3), pointer :: string_nl(:)
        character(len=1)       :: l1(0:20)

        ! Initial values
        data l1 /"s","p","d","f","g","h","i","k","l","m","n", &
            "o","q","r","t","u","v","w","x","y","z"/
        counter = 0
        allocate(occupation(csf_set_ls%shell_count))
        allocate(string_ls(csf_set_ls%shell_count))
        allocate(string_xls(csf_set_ls%shell_count))
        do j = csf_set_ls%closed_shell_count+1, csf_set_ls%shell_count
            if (csf_set_ls%csf(csf_index)%occupation(j) > 0) then
                counter = counter +1;        occupation(counter) = j
                call format_ls_csf_label(csf_index,j,ls,xls)
                string_ls(counter)  = ls;    string_xls(counter) = xls
            end if
        end do
        allocate(string_nl(counter))
        do i = 1,counter
            write(string,"(I2)")csf_set_ls%shell(occupation(i))%n
            string_nl(i)(1:2) = string(1:2)
            string_nl(i)(3:3) = l1(csf_set_ls%shell(occupation(i))%l)
        end do
        if(counter > 18) then
            print*, counter
            stop "pack_ls_csf_label: A"
        end if
        do i = 1,counter
            couple(i)(1:3) = string_ls(i)(1:3)
            if(i > 1) then
                couple(i+counter-1)(1:3) = string_xls(i)(1:3)
            end if
            q(i) = csf_set_ls%csf(csf_index)%occupation(occupation(i))
            elc(i)(1:3) = string_nl(i)(1:3)
        end do
        call pack_ls_subshell_label(counter,elc,q,couple,string_csf_pack)
        string_csf = string_csf_pack
        deallocate(string_nl)
        deallocate(occupation)
        deallocate(string_ls)
        deallocate(string_xls)
    end subroutine pack_ls_csf_label

    !     Print all information about the single CSF scheme csf_set in
    !     a nead format on stream.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: Dec 2015
    subroutine write_jj_csf(stream,csf_index,csf_number_local)
        ! Module dependencies
        use orbital_csf_basis,           only: orbital_count, orbital_principal_n, orbital_kappa
        use csf_pair_angular_state,             only: core_orbital_count
        use csf_quantum_queries, only: csf_total_j_code, csf_occupation, csf_subshell_quantum_code
        use half_integer_format, only: format_half_integer
        implicit none (type, external)

        ! Arguments and local declarations
        save  ! Preserve the legacy serial controller workspace.
        integer, intent(in) :: stream, csf_index, csf_number_local
        integer                                 :: i, orbital_index, counter, text_length
        integer, pointer          :: occupation(:)
        character(len=4), pointer :: string_j(:), string_x(:)
        character(len=4)                        :: string_cnum
        character(len=2)                        :: string
        character(len=4), pointer ::  string_nl(:)
        character(len=1)       :: l1(0:20)

        ! Initial values
        data l1 /"s","p","d","f","g","h","i","k","l","m","n", &
            "o","q","r","t","u","v","w","x","y","z"/
        counter = 0
        allocate(occupation(orbital_count))
        allocate(string_j(orbital_count))
        allocate(string_x(orbital_count))
        do orbital_index = core_orbital_count+1, orbital_count
            if (csf_occupation(orbital_index,csf_index) > 0) then
                counter = counter +1
                occupation(counter) = orbital_index
                if (csf_subshell_quantum_code(3,orbital_index,csf_index)-1 == 0) then
                    string_j(counter) = "    "
                else
                    call format_half_integer(csf_subshell_quantum_code(3,orbital_index,csf_index)-1,string_cnum,text_length)
                    string_j(counter) = string_cnum(1:text_length)
                end if
                if (intermediate_two_j(orbital_index,csf_index) == 0) then
                    string_x(counter) = "     "
                else
                    call format_half_integer(intermediate_two_j(orbital_index,csf_index),string_cnum,text_length)
                    string_x(counter) =  string_cnum(1:text_length)
                end if
            end if
        end do
        allocate(string_nl(counter))
        do i = 1,counter
            write(string,"(I2)")orbital_principal_n(occupation(i))
            string_nl(i)(1:2) = string(1:2)
            orbital_index = ((iabs(orbital_kappa(occupation(i)))*2)-1+    &
                orbital_kappa(occupation(i))/iabs(orbital_kappa(occupation(i))))/2
            string_nl(i)(3:3) = l1(orbital_index)
            if(orbital_kappa(occupation(i)) <= 0) then
                string_nl(i)(4:4) = " "
            else
                string_nl(i)(4:4) = "-"
            end if
        end do
        if (stream == -1) then
            write(*,'(4x,i9,")",100(a4,"(",i2,")",2x))') csf_number_local,(string_nl(orbital_index),                   &
                csf_occupation(occupation(orbital_index),csf_index),orbital_index=1,counter)
            write(*,'(18x,a4,100(6x,a4))')(string_j(orbital_index),orbital_index=1,counter)
            call format_half_integer(csf_total_j_code(csf_index)-1,string_cnum,text_length)
            write(*,'(31x,a4,100(6x,a5))')(string_x(orbital_index),orbital_index=2,counter-1),string_cnum(1:text_length)
        end if
        deallocate(occupation)
        deallocate(string_j)
        deallocate(string_x)
    end subroutine write_jj_csf

    !     Print all information about the single CSF scheme csf_set in
    !     a nead format on stream.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: May 2011
    subroutine write_ls_csf(stream,csf_set_ls,csf_index)
        ! Module dependencies
        use half_integer_format, only: format_half_integer
        implicit none (type, external)

        ! Arguments and local declarations
        save  ! Preserve the legacy serial controller workspace.
        integer, intent(in)            :: stream, csf_index
        type(csf_basis_ls), intent(in) :: csf_set_ls
        integer                                 :: i, j, counter, text_length
        integer, pointer          :: occupation(:)
        character(len=4)                        :: ls, xls
        character(len=4),pointer  :: string_ls(:), string_xls(:)
        character(len=4)                        :: string_cnum
        character(len=2)                        :: string
        character(len=3), pointer :: string_nl(:)
        character(len=1)       :: l1(0:20)

        ! Initial values
        data l1 /"s","p","d","f","g","h","i","k","l","m","n", &
            "o","q","r","t","u","v","w","x","y","z"/
        counter = 0
        allocate(occupation(csf_set_ls%shell_count))
        allocate(string_ls(csf_set_ls%shell_count))
        allocate(string_xls(csf_set_ls%shell_count))
        do j = csf_set_ls%closed_shell_count+1, csf_set_ls%shell_count
            if (csf_set_ls%csf(csf_index)%occupation(j) > 0) then
                counter = counter +1;        occupation(counter) = j
                call format_ls_csf(csf_index,j,ls,xls)
                string_ls(counter)  = ls;    string_xls(counter) = xls
            end if
        end do
        allocate(string_nl(counter))
        do i = 1,counter
            write(string,"(I2)")csf_set_ls%shell(occupation(i))%n
            string_nl(i)(1:2) = string(1:2)
            string_nl(i)(3:3) = l1(csf_set_ls%shell(occupation(i))%l)
        end do
        if (stream == -1) then
            write(*,'(i10,")",4x,100(a3,"(",i2,")",2x))')csf_index,(string_nl(j),                                 &
                csf_set_ls%csf(csf_index)%occupation(occupation(j)),j=1,counter)
            call format_half_integer(csf_set_ls%csf(csf_index)%two_j,string_cnum,text_length)
            write(*,'(19x,a4,100(5x,a4))')(string_ls(j),j=1,counter),(string_xls(j),j=2,counter),&
                string_cnum(1:text_length)
        else
            write(stream,'(8(1x,a3,"(",i2,")"))')(string_nl(j),                         &
                csf_set_ls%csf(csf_index)%occupation(occupation(j)),j=1,counter)
            call format_half_integer(csf_set_ls%csf(csf_index)%two_j,string_cnum,text_length)
            write(stream,'(1x,15(a4))')(string_ls(j),j=1,counter),(string_xls(j),j=2,counter)
        end if
        deallocate(string_nl)
        deallocate(occupation)
        deallocate(string_ls)
        deallocate(string_xls)
    end subroutine write_ls_csf

    ! Print all information about the CSF scheme csf_set in a nead
    ! format on stream.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: May 2011
    subroutine write_ls_csf_list(stream,csf_set_ls,before_j)
        implicit none (type, external)

        ! Arguments and local declarations
        save  ! Preserve the legacy serial controller workspace.
        integer, intent(in)            :: stream
        type(csf_basis_ls), intent(in) :: csf_set_ls
        integer, intent(in)            :: before_j
        integer                                  :: i, j
        character (len=2)                        :: string
        character (len=3), pointer :: string_nl(:)
        character (len=1)       :: l1(0:20)

        ! Initial values
        data l1 /"s","p","d","f","g","h","i","k","l","m","n", &
            "o","q","r","t","u","v","w","x","y","z"/
        if(before_j == 0) then
            write(stream,*) " "
            allocate(string_nl(csf_set_ls%closed_shell_count))
            do i = 1,csf_set_ls%closed_shell_count
                write(string,"(I2)")csf_set_ls%shell(i)%n
                string_nl(i)(1:2) = string(1:2)
                string_nl(i)(3:3) = l1(csf_set_ls%shell(i)%l)
            end do
            write(stream,'(18(1x,a3))')(string_nl(j),j=1,csf_set_ls%closed_shell_count)
            deallocate(string_nl)
        end if
        ! output to *.lsj.c
        do i = 1,csf_set_ls%ls_csf_count
            call write_ls_csf(stream,csf_set_ls,i)
        end do
        write(stream,*) "*"
    end subroutine write_ls_csf_list

    !     This subroutine fills up the variable asf_set_LS%csf_set_LS
    !     with data generated using the one from asf_set%csf_set.
    !
    !     This subroutine contains the following internal subroutines:
    !          * subroutine            setLS_action
    !          * subroutine            setLS_add_quantum_numbers
    !          * recursive subroutine  setLS_job_count
    !          * function              setLS_equivalent_csfs
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: Dec 2015
    subroutine build_ls_basis(retained_csf_mask,first_csf_index,last_csf_index)
        ! Module dependencies
        use atomic_calculation_state,           only: electron_count
        use orbital_csf_basis,           only: orbital_count, orbital_principal_n, orbital_kappa, csf_count
        use csf_pair_angular_state,             only: core_orbital_count
        use csf_quantum_queries, only: csf_total_j_code, csf_occupation
        implicit none (type, external)

        ! Arguments and local declarations
        save  ! Preserve the legacy serial controller workspace.
        integer, intent(in)  :: retained_csf_mask(:)
        integer, intent(in)  :: first_csf_index, last_csf_index
        type(nonrelativistic_subshell_label), pointer :: shell_temp(:)
        integer, pointer  :: nllsval(:), nljjval(:)
        type(lsj_list)                  :: nonequiv_csfs_jj
        integer :: isubc, isubc2, icsf_jj, icsf_jj2, icsf_jj_real, icsf_ls
        integer :: n, action_type
        logical :: new_one, found_parent_minus, found_parent_plus
        integer, pointer :: all_occupation(:)
        integer, pointer :: li(:), si(:), l_i(:), s_i(:), w(:), q(:)
        integer                        :: j
        asf_set_ls%csf_set_ls%number_of_electrons = electron_count
        ! 1. define  nl, parent
        allocate(shell_temp(orbital_count));  allocate(nllsval(orbital_count));  allocate(nljjval(orbital_count))
        nljjval = 0;  nllsval = 0
        asf_set_ls%csf_set_ls%shell_count = 0
        do isubc = 1, orbital_count, 1
            nljjval(isubc) = orbital_principal_n(isubc)*100 +                              &
                ((iabs(orbital_kappa(isubc))*2)-1+orbital_kappa(isubc)/iabs(orbital_kappa(isubc)))/2
            new_one = .true.
            do isubc2 = 1, isubc-1, 1
                if(nljjval(isubc) == nljjval(isubc2)) new_one = .false.
            end do
            if(new_one) then
                asf_set_ls%csf_set_ls%shell_count = asf_set_ls%csf_set_ls%shell_count + 1
                shell_temp(asf_set_ls%csf_set_ls%shell_count)%n = orbital_principal_n(isubc)
                shell_temp(asf_set_ls%csf_set_ls%shell_count)%l=                &
                    ((iabs(orbital_kappa(isubc))*2)-1+orbital_kappa(isubc)/iabs(orbital_kappa(isubc)))/2
                nllsval(asf_set_ls%csf_set_ls%shell_count) = nljjval(isubc)
            end if
        end do
        allocate(asf_set_ls%csf_set_ls%shell(asf_set_ls%csf_set_ls%shell_count))
        allocate(asf_set_ls%csf_set_ls%parent(asf_set_ls%csf_set_ls%shell_count))
        do isubc=1, asf_set_ls%csf_set_ls%shell_count, 1
            asf_set_ls%csf_set_ls%shell(isubc) = shell_temp(isubc)
        end do
        deallocate(shell_temp)
        ! begin find parent
        do isubc = 1, asf_set_ls%csf_set_ls%shell_count, 1
            found_parent_minus = .false.;  found_parent_plus = .false.
            do isubc2 = 1, orbital_count, 1
                if(nllsval(isubc) == nljjval(isubc2)) then
                    if(orbital_kappa(isubc2) > 0) then
                        found_parent_minus = .true.
                        asf_set_ls%csf_set_ls%parent(isubc)%parent_minus = isubc2
                    else
                        found_parent_plus = .true.
                        asf_set_ls%csf_set_ls%parent(isubc)%parent_plus  = isubc2
                    end if
                end if
            end do
            if(.not.found_parent_plus)  &
                asf_set_ls%csf_set_ls%parent(isubc)%parent_plus = 0
            if(.not.found_parent_minus)  &
                asf_set_ls%csf_set_ls%parent(isubc)%parent_minus = 0
        end do
        ! end find parent  --------------------
        !
        ! 2. define the number of "core" shells
        !    (the LS shell is supposed to be "core" if:
        !    1. l=0 and corresponding jj subshell is "core"
        !    2. l<>0 l+ and l - "core" subshells            )
        asf_set_ls%csf_set_ls%closed_shell_count=0
        do isubc=1, asf_set_ls%csf_set_ls%shell_count, 1
            if(asf_set_ls%csf_set_ls%parent(isubc)%parent_minus <= core_orbital_count &
                .and.                                                       &
                asf_set_ls%csf_set_ls%parent(isubc)%parent_plus <= core_orbital_count) &
                asf_set_ls%csf_set_ls%closed_shell_count = asf_set_ls%csf_set_ls%closed_shell_count+1
        end do
        ! 3. form the list of "nonequivalent" csfs_jj
        !    (i.e. csfs_jj different in J,parity, or
        !     some l's occupation numbers Ni = N_(i+) + N_(i-))
        allocate(nonequiv_csfs_jj%items(csf_count))
        nonequiv_csfs_jj%list_size = 0
        do icsf_jj = first_csf_index, last_csf_index
            new_one = .true.
            if(retained_csf_mask(icsf_jj) == 1) then
                do icsf_jj2 = first_csf_index , icsf_jj - 1
                    if(retained_csf_mask(icsf_jj2) == 1) then
                        if(are_ls_csfs_equivalent(icsf_jj,icsf_jj2)) then
                            new_one = .false.
                            exit
                        end if
                    end if
                end do
                if(new_one) then
                    nonequiv_csfs_jj%list_size = nonequiv_csfs_jj%list_size + 1
                    nonequiv_csfs_jj%items(nonequiv_csfs_jj%list_size) = icsf_jj
                end if
            end if
        end do
        ! 4. for each nonequivalent csf_jj find all the csfs_LS
        !
        !     To avoid the dependency on the number of subshells
        !     the recursive subroutine is used
        allocate(li(asf_set_ls%csf_set_ls%shell_count))
        allocate(l_i(asf_set_ls%csf_set_ls%shell_count))
        allocate(si(asf_set_ls%csf_set_ls%shell_count))
        allocate(s_i(asf_set_ls%csf_set_ls%shell_count))
        allocate(q(asf_set_ls%csf_set_ls%shell_count))
        allocate(w(asf_set_ls%csf_set_ls%shell_count))
        allocate(all_occupation(asf_set_ls%csf_set_ls%shell_count))
        ! 4.1 - find the number of csfs_LS
        asf_set_ls%csf_set_ls%ls_csf_count = 0
        do icsf_jj=1, nonequiv_csfs_jj%list_size, 1
            ! set the variable for conviniency ...
            icsf_jj_real = nonequiv_csfs_jj%items(icsf_jj)
            j = csf_total_j_code(icsf_jj_real)-1
            ! define the occupation numbers
            do isubc = 1, asf_set_ls%csf_set_ls%shell_count, 1
                all_occupation(isubc) = 0
                do isubc2 = 1, orbital_count, 1
                    if(nllsval(isubc) == nljjval(isubc2))       &
                        all_occupation(isubc)=all_occupation(isubc)+csf_occupation(isubc2,icsf_jj_real)
                end do  !isubc2
            end do  !isubc
            n = 0;    action_type = 1
            call count_ls_basis_jobs(1, n)
            asf_set_ls%csf_set_ls%ls_csf_count=asf_set_ls%csf_set_ls%ls_csf_count + n
        end do
        ! 4.2 - fill up the csf_LS arrays with the corresponding quantum numbers
        ! 4.2.1  allocate the arrays
        allocate(asf_set_ls%csf_set_ls%csf(asf_set_ls%csf_set_ls%ls_csf_count))
        do n = 1, asf_set_ls%csf_set_ls%ls_csf_count, 1
            allocate(asf_set_ls%csf_set_ls%csf(n)%occupation(asf_set_ls%csf_set_ls%shell_count))
            allocate(asf_set_ls%csf_set_ls%csf(n)%seniority(asf_set_ls%csf_set_ls%shell_count))
            allocate(asf_set_ls%csf_set_ls%csf(n)%w(asf_set_ls%csf_set_ls%shell_count))
            allocate(asf_set_ls%csf_set_ls%csf(n)%shell_two_l(asf_set_ls%csf_set_ls%shell_count))
            allocate(asf_set_ls%csf_set_ls%csf(n)%shell_two_s(asf_set_ls%csf_set_ls%shell_count))
            allocate(asf_set_ls%csf_set_ls%csf(n)%intermediate_two_l(asf_set_ls%csf_set_ls%shell_count))
            allocate(asf_set_ls%csf_set_ls%csf(n)%intermediate_two_s(asf_set_ls%csf_set_ls%shell_count))
        end do
        ! 4.2.2 - fill them with quantum numbers
        icsf_ls = 0
        do icsf_jj=1, nonequiv_csfs_jj%list_size, 1
            ! set the variable for conveniency ...
            icsf_jj_real = nonequiv_csfs_jj%items(icsf_jj)
            j = csf_total_j_code(icsf_jj_real)-1
            do isubc = 1, asf_set_ls%csf_set_ls%shell_count, 1
                all_occupation(isubc)=0
                do isubc2 = 1, orbital_count, 1
                    if(nllsval(isubc) == nljjval(isubc2))            &
                        all_occupation(isubc)=all_occupation(isubc)+csf_occupation(isubc2,icsf_jj_real)
                end do  !isubc2
            end do  !isubc
            action_type = 2
            call count_ls_basis_jobs(1, n)
        end do
        deallocate(nonequiv_csfs_jj%items)
        deallocate(nllsval);  deallocate(nljjval)
        deallocate(q);   deallocate(w);    deallocate(all_occupation)
        deallocate(li);  deallocate(l_i);  deallocate(si)
        deallocate(s_i)
    contains

        ! The subroutine defines the "action" of subroutine
        ! setLS_job_count: whether it counts the number
        ! of csfs_LS (asf_set_LS%csf_set_LS%novcsf) or fills the
        ! arrays of wave functions in LS coupling with
        ! asf_set_LS%csf_set_LS%csf(...) with the corresponding
        ! quantum nubers.
        !
        ! Written by G. Gaigalas,
        ! NIST                                     last update: May 2011
        subroutine apply_ls_basis_action(basis_action,generated_csf_count)
            implicit none (type, external)

            ! Local declarations
            save  ! Preserve the legacy serial controller workspace.

            ! Arguments
            integer :: basis_action, generated_csf_count
            if (basis_action == 1) then
                generated_csf_count = generated_csf_count + 1
            else
                call add_ls_quantum_numbers()
            end if
        end subroutine apply_ls_basis_action

        !     The subroutine adds quantum numbers stored in temprorary arrays
        !     Li, Si, L_i, S_i, w, Q to the corresponding arrays
        !     of asf_set_LS%csf_set_LS%csf() (i.e. to the arrays of
        !     corresponding quantum numbers of the wave function in LS
        !     coupling).
        !
        ! Written by G. Gaigalas,
        ! NIST                                     last update: May 2011
        subroutine add_ls_quantum_numbers()
            ! Module dependencies
            use csf_quantum_queries, only: csf_parity_code
            implicit none (type, external)

            ! Local declarations
            save  ! Preserve the legacy serial controller workspace.
            integer :: shell_index
            icsf_ls = icsf_ls + 1
            if(icsf_ls  >  asf_set_ls%csf_set_ls%ls_csf_count) then
                stop "add_ls_quantum_numbers(): program stop A."
            end if
            asf_set_ls%csf_set_ls%csf(icsf_ls)%two_j = j
            if (csf_parity_code(icsf_jj_real) == 1) then
                asf_set_ls%csf_set_ls%csf(icsf_ls)%parity = "+"
            else if (csf_parity_code(icsf_jj_real) == -1) then
                asf_set_ls%csf_set_ls%csf(icsf_ls)%parity = "-"
            end if
            do shell_index = 1, asf_set_ls%csf_set_ls%shell_count, 1
                asf_set_ls%csf_set_ls%csf(icsf_ls)%occupation(shell_index) =        &
                    all_occupation(shell_index)
                asf_set_ls%csf_set_ls%csf(icsf_ls)%shell_two_l(shell_index)  = li(shell_index)
                asf_set_ls%csf_set_ls%csf(icsf_ls)%shell_two_s(shell_index)  = si(shell_index)
                asf_set_ls%csf_set_ls%csf(icsf_ls)%intermediate_two_l(shell_index) = l_i(shell_index)
                asf_set_ls%csf_set_ls%csf(icsf_ls)%intermediate_two_s(shell_index) = s_i(shell_index)
                asf_set_ls%csf_set_ls%csf(icsf_ls)%w(shell_index)       = w(shell_index)
                asf_set_ls%csf_set_ls%csf(icsf_ls)%seniority(shell_index) =         &
                    2*asf_set_ls%csf_set_ls%shell(shell_index)%l + 1 - q(shell_index)
            end do
        end subroutine add_ls_quantum_numbers

        !     Recursive subroutine for the calculation of the
        !     number of csfs_LS and corresponding quantum numbers.
        !
        ! Written by G. Gaigalas,
        ! NIST                                     last update: Dec 2015
        recursive subroutine count_ls_basis_jobs(isubc, rez)
            ! Module dependencies
            use ls_jj_subshell_terms, only: subshell_term_ls
            use angular_triangle_and_parity, only: satisfies_angular_triangle_and_parity
            implicit none (type, external)

            ! Arguments
            integer :: isubc, rez

            ! Local declarations
            type(subshell_term_ls) :: ls_terms(120)
            integer :: iterm
            ! GG 2015.12.02 NIST
            integer :: number, numbergg
            integer :: templmax,templmin,tempsmax,tempsmin,temps,templ
            integer :: l_shell, n
            if(isubc>(asf_set_ls%csf_set_ls%shell_count) .or. isubc<1) then
                print *, "isubc = ", isubc
                stop "count_ls_basis_jobs(): program stop A."
            end if
            if(isubc<=asf_set_ls%csf_set_ls%shell_count) then
                if(all_occupation(isubc)==0) then
                    if(isubc>1) then
                        li(isubc)  = 0;          l_i(isubc) = l_i(isubc-1)
                        si(isubc)  = 0;          s_i(isubc) = s_i(isubc-1)
                        if(isubc < asf_set_ls%csf_set_ls%shell_count) then
                            call count_ls_basis_jobs(isubc + 1, rez)
                        else
                            if(satisfies_angular_triangle_and_parity(s_i(isubc),l_i(isubc),j)==1)          &
                                call apply_ls_basis_action(action_type, rez) !rez=rez+1
                        end if
                    else
                        li(isubc)  = 0;          l_i(isubc) = 0
                        si(isubc)  = 0;          s_i(isubc) = 0
                        !  Modification for configurations with empty 1s subshell
                        call count_ls_basis_jobs(isubc + 1, rez)
                        ! CYC-2023 end
                    end if
                else
                    n = all_occupation(isubc)
                    l_shell = asf_set_ls%csf_set_ls%shell(isubc)%l
                    call enumerate_ls_subshell_terms(l_shell,n,ls_terms,number)
                    do iterm=1, number, 1
                        ! GG 2015.12.02 NIST
                        n = all_occupation(isubc)
                        l_shell = asf_set_ls%csf_set_ls%shell(isubc)%l
                        call enumerate_ls_subshell_terms(l_shell,n,ls_terms,numbergg)
                        ! GG 2015.12.02 NIST
                        li(isubc)  = ls_terms(iterm)%subshell_two_l
                        si(isubc)  = ls_terms(iterm)%subshell_two_s
                        w(isubc)   = ls_terms(iterm)%w
                        q(isubc)   = ls_terms(iterm)%subshell_two_quasispin
                        if(isubc==1) then
                            l_i(isubc) = ls_terms(iterm)%subshell_two_l
                            s_i(isubc) = ls_terms(iterm)%subshell_two_s
                            if(asf_set_ls%csf_set_ls%shell_count>1) then
                                call count_ls_basis_jobs(isubc + 1, rez)
                            else
                                if(satisfies_angular_triangle_and_parity(s_i(isubc),l_i(isubc),j)==1)          &
                                    call apply_ls_basis_action(action_type, rez) !rez=rez+1
                            end if
                        else
                            templmax=l_i(isubc-1)+li(isubc)
                            templmin=abs(l_i(isubc-1)-li(isubc))
                            tempsmax=s_i(isubc-1)+si(isubc)
                            tempsmin=abs(s_i(isubc-1)-si(isubc))
                            do templ = templmin, templmax, 2
                                l_i(isubc) = templ
                                do temps = tempsmin, tempsmax, 2
                                    s_i(isubc) = temps
                                    if(isubc<asf_set_ls%csf_set_ls%shell_count) then
                                        call count_ls_basis_jobs(isubc+1,rez)
                                    else
                                        if(satisfies_angular_triangle_and_parity(s_i(isubc),l_i(isubc),j)==1)        &
                                            call apply_ls_basis_action(action_type, rez) ! rez=rez+1
                                    end if
                                end do   ! tempS
                            end do      ! tempL
                        end if
                    end do            ! iterm
                end if
            end if
        end subroutine count_ls_basis_jobs

        !     This subroutine defines the "equivalency" of two csfs_jj
        !     in the sence of generation of csfs_LS
        !     number of csfs_LS and corresponding quantum numbers.
        !
        ! Written by G. Gaigalas,
        ! NIST                                     last update: May 2011
        function are_ls_csfs_equivalent(first_csf_index,second_csf_index)         result(equivalent)
            ! Module dependencies
            use orbital_csf_basis,           only: orbital_count
            use csf_quantum_queries, only: csf_parity_code
            implicit none (type, external)

            ! Arguments
            integer :: first_csf_index, second_csf_index
            logical :: equivalent

            ! Local declarations
            integer :: isubc_ls, orbital_index, nls1, nls2
            equivalent= .true.
            if(csf_total_j_code(first_csf_index) /= csf_total_j_code(second_csf_index)) then
                equivalent = .false.
            else if (csf_parity_code(first_csf_index) /= csf_parity_code(second_csf_index)) then
                equivalent = .false.
            else
                do isubc_ls = asf_set_ls%csf_set_ls%closed_shell_count + 1,  &
                    asf_set_ls%csf_set_ls%shell_count, 1
                    nls1 = 0;   nls2 = 0
                    do orbital_index = 1, orbital_count, 1
                        if(nllsval(isubc_ls) == nljjval(orbital_index)) then
                            nls1 = nls1 + csf_occupation(orbital_index,first_csf_index)
                            nls2 = nls2 + csf_occupation(orbital_index,second_csf_index)
                        end if
                    end do
                    if(nls1/=nls2) then
                        equivalent=.false.
                        return
                    end if
                end do
            end if
        end function are_ls_csfs_equivalent
    end subroutine build_ls_basis

    !     This procedure return the value of the transformation matrix
    !     from jj- to LS-coupling scheme in the case of any number of
    !     open shells.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: May 2011
    recursive function ls_jj_csf_overlap(jj_csf_index,ls_csf_index)          result(coefficient_value)
        ! Module dependencies
        use orbital_csf_basis,           only: orbital_kappa
        use numeric_constants,          only: zero
        use csf_quantum_queries, only: csf_occupation, csf_subshell_quantum_code
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: jj_csf_index,ls_csf_index
        real(real64)        :: coefficient_value

        ! Local declarations
        integer :: total_number, shell_number, number_minus, number_plus, &
            jj_minus, n_minus, q_minus, j_1_minus,                 &
            jj_plus, n_plus, q_plus, j_1_plus, j_1,                &
            l_shell, n_ls, w_1, q_1, l_1, s_1
        shell_number=asf_set_ls%csf_set_ls%closed_shell_count+1
        number_minus=asf_set_ls%csf_set_ls%parent(shell_number)%parent_minus
        number_plus =asf_set_ls%csf_set_ls%parent(shell_number)%parent_plus
        if (number_minus+1 /= number_plus .and.  &
            number_minus*number_plus  /= 0) then
            stop "ls_jj_csf_overlap(): program stop A."
        end if
        if (number_minus == 0) then
            jj_minus = iabs(orbital_kappa(number_plus))*2 - 3
            n_minus  = 0;    j_1_minus = 0;   q_minus = (jj_minus + 1)/2
        else
            jj_minus  = iabs(orbital_kappa(number_minus))*2 - 1
            n_minus   = csf_occupation(number_minus,jj_csf_index)
            q_minus   = (jj_minus +1)/2 - csf_subshell_quantum_code(1,number_minus,jj_csf_index)
            j_1_minus = csf_subshell_quantum_code(3,number_minus,jj_csf_index)-1
            j_1       = intermediate_two_j(number_minus,jj_csf_index)
        end if
        if (number_plus == 0) then
            jj_plus = iabs(orbital_kappa(number_minus))*2 + 1
            n_plus  = 0;    j_1_plus = 0;   q_plus = (jj_plus + 1)/2
        else
            jj_plus  = iabs(orbital_kappa(number_plus))*2 - 1
            n_plus   = csf_occupation(number_plus,jj_csf_index)
            q_plus   = (jj_plus + 1)/2 - csf_subshell_quantum_code(1,number_plus,jj_csf_index)
            j_1_plus = csf_subshell_quantum_code(3,number_plus,jj_csf_index)-1
            j_1       = intermediate_two_j(number_plus,jj_csf_index)
        end if
        l_shell=asf_set_ls%csf_set_ls%shell(shell_number)%l
        n_ls   =asf_set_ls%csf_set_ls%csf(ls_csf_index)%occupation(shell_number)
        w_1    =asf_set_ls%csf_set_ls%csf(ls_csf_index)%w(shell_number)
        q_1    =2*l_shell+1-                                              &
            asf_set_ls%csf_set_ls%csf(ls_csf_index)%seniority(shell_number)
        l_1     = asf_set_ls%csf_set_ls%csf(ls_csf_index)%shell_two_l(shell_number)
        s_1     = asf_set_ls%csf_set_ls%csf(ls_csf_index)%shell_two_s(shell_number)
        coefficient_value = ls_jj_recoupling_coefficient(l_shell,n_ls,w_1,q_1,l_1,s_1,j_1,                   &
            jj_minus,n_minus,q_minus,j_1_minus,   &
            jj_plus,q_plus,j_1_plus)
        if (abs(coefficient_value) > coefficient_threshold) then
            total_number =                                                 &
                asf_set_ls%csf_set_ls%shell_count - asf_set_ls%csf_set_ls%closed_shell_count
            if (total_number == 1) then
                coefficient_value = coefficient_value
            else if (total_number >= 2) then
                do shell_number = asf_set_ls%csf_set_ls%closed_shell_count + 2,         &
                    asf_set_ls%csf_set_ls%shell_count
                    if (abs(coefficient_value) > coefficient_threshold) then
                        coefficient_value = coefficient_value * ls_jj_shell_recoupling_factor(shell_number,jj_csf_index, &
                            & ls_csf_index)
                    end if
                end do
            else
                coefficient_value = zero
            end if
        else
            coefficient_value = zero
        end if
    end function ls_jj_csf_overlap

    !     Return the value of main part of the transformation matrix from
    !     jj- to LS-coupling scheme in the case of any number of open
    !     shells.
    !
    ! Written by G. Gaigalas,
    ! NIST                                     last update: May 2011
    ! Modified by G. Gaigalas. Jan 2025
    recursive function ls_jj_shell_recoupling_factor(shell_number,jj_csf_index,ls_csf_index)   result(coefficient_value)
        ! Module dependencies
        use orbital_csf_basis,           only: orbital_kappa
        use numeric_constants,          only: zero, one
        use wigner_6j_triangle_selection, only: allows_wigner_6j_triangles
        use csf_quantum_queries, only: csf_occupation, csf_subshell_quantum_code
        use wigner_9j, only: wigner_9j_value
        use wigner_6j, only: wigner_6j_value
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: shell_number,jj_csf_index,ls_csf_index
        real(real64)        :: coefficient_value

        ! Local declarations
        real(real64) :: wa_sum, rac6, rac9
        integer      :: delta_j, number_minus, number_plus,   &
            ! GG begin  2025
            number_plus_1,                        &
            ! GG end    2025
            jj_minus, n_minus, q_minus,           &
            jj_plus, n_plus, q_plus,              &
            j_i_min, j_i_max, j_i_minus,          &
            j_i_plus, j_i, j_1_i, jp_1_i, j_1_i1, &
            l_shell, n_ls, w_i, q_i, l_i, s_i,    &
            l_1_i, s_1_i, l_1_i1, s_1_i1
        coefficient_value     = zero;      wa_sum = zero
        number_minus = asf_set_ls%csf_set_ls%parent(shell_number)%parent_minus
        number_plus  = asf_set_ls%csf_set_ls%parent(shell_number)%parent_plus
        if (number_minus+1 /= number_plus .and.  &
            number_minus*number_plus  /= 0) then
            stop "tranLSjjmp(): program stop A."
        end if
        number_plus_1=asf_set_ls%csf_set_ls%parent(shell_number-1)%parent_plus
        if (number_plus_1 == 0) then
            number_plus_1 =                                            &
                asf_set_ls%csf_set_ls%parent(shell_number-1)%parent_minus
        end if
        if (number_minus == 0) then
            jj_minus  = iabs(orbital_kappa(number_plus))*2 - 1
            n_minus   = 0;    j_i_minus = 0;    q_minus = (jj_minus + 1)/2
            jp_1_i    = intermediate_two_j(number_plus_1,jj_csf_index)
        else
            jj_minus  = iabs(orbital_kappa(number_minus))*2 - 1
            n_minus   = csf_occupation(number_minus,jj_csf_index)
            q_minus   = (jj_minus + 1)/2 - csf_subshell_quantum_code(1,number_minus,jj_csf_index)
            j_i_minus = csf_subshell_quantum_code(3,number_minus,jj_csf_index)-1
            jp_1_i    = intermediate_two_j(number_minus,jj_csf_index)
        end if
        if (number_plus == 0) then
            jj_plus  = iabs(orbital_kappa(number_minus))*2 - 1
            q_plus   = (jj_plus + 1)/2
            n_plus   = 0;         j_i_plus = 0;         j_1_i   = jp_1_i
        else
            jj_plus  = iabs(orbital_kappa(number_plus))*2 - 1
            n_plus   = csf_occupation(number_plus,jj_csf_index)
            q_plus   = (jj_plus + 1)/2 - csf_subshell_quantum_code(1,number_plus,jj_csf_index)
            j_i_plus = csf_subshell_quantum_code(3,number_plus,jj_csf_index)-1
            j_1_i     = intermediate_two_j(number_plus,jj_csf_index)
        end if
        j_i_min = iabs(j_i_minus - j_i_plus);   j_i_max = j_i_minus + j_i_plus
        l_shell = asf_set_ls%csf_set_ls%shell(shell_number)%l
        n_ls    = asf_set_ls%csf_set_ls%csf(ls_csf_index)%occupation(shell_number)
        w_i     = asf_set_ls%csf_set_ls%csf(ls_csf_index)%w(shell_number)
        q_i     = 2 * l_shell + 1 -                                          &
            asf_set_ls%csf_set_ls%csf(ls_csf_index)%seniority(shell_number)
        ! GG begin  2025
        if (n_ls == n_minus + n_plus .and. n_ls == 0) then
            coefficient_value = one
        else if (n_ls == n_minus + n_plus .and. n_ls /= 0) then
            ! GG end    2025
            l_i = asf_set_ls%csf_set_ls%csf(ls_csf_index)%shell_two_l(shell_number)
            s_i = asf_set_ls%csf_set_ls%csf(ls_csf_index)%shell_two_s(shell_number)
            l_1_i = asf_set_ls%csf_set_ls%csf(ls_csf_index)%intermediate_two_l(shell_number)
            s_1_i = asf_set_ls%csf_set_ls%csf(ls_csf_index)%intermediate_two_s(shell_number)
            if (shell_number == 2) then
                ! GG begin  2025
                j_1_i1 = intermediate_two_j(number_plus_1,jj_csf_index)
                ! GG end    2025
                l_1_i1 =asf_set_ls%csf_set_ls%csf(ls_csf_index)%shell_two_l(shell_number-1)
                s_1_i1 =asf_set_ls%csf_set_ls%csf(ls_csf_index)%shell_two_s(shell_number-1)
            else
                j_1_i1 =intermediate_two_j(number_plus_1,jj_csf_index)
                l_1_i1 =asf_set_ls%csf_set_ls%csf(ls_csf_index)%intermediate_two_l(shell_number-1)
                s_1_i1 =asf_set_ls%csf_set_ls%csf(ls_csf_index)%intermediate_two_s(shell_number-1)
            end if
            do j_i = j_i_min, j_i_max, 2
                delta_j = allows_wigner_6j_triangles(j_i_minus,j_i_plus,j_i,j_1_i,j_1_i1,jp_1_i)
                if (delta_j /= 0) then
                    call wigner_9j_value(l_1_i1,s_1_i1,j_1_i1,l_i,s_i,j_i,l_1_i,s_1_i,j_1_i,  &
                        1,delta_j,rac9)
                    if (delta_j /= 0) then
                        call wigner_6j_value(j_i_minus,j_i_plus,j_i,j_1_i,j_1_i1,jp_1_i,0,rac6)
                        call wigner_9j_value(l_1_i1,s_1_i1,j_1_i1,l_i,s_i,j_i,l_1_i,s_1_i,j_1_i,&
                            0,delta_j,rac9)
                        wa_sum = wa_sum + (j_i + one) * rac6 * rac9 *                &
                            ls_jj_recoupling_coefficient(l_shell,n_ls,w_i,q_i,l_i,s_i,j_i,              &
                            jj_minus,n_minus,q_minus,j_i_minus,    &
                            jj_plus,q_plus,j_i_plus)
                    end if
                end if
            end do
            coefficient_value = wa_sum * sqrt((jp_1_i+one)*(j_1_i1+one)*(l_1_i+one)*(s_1_i+one))
            if (mod(j_i_minus+j_i_plus+j_1_i1+j_1_i,4) /= 0) coefficient_value = - coefficient_value
        else
            coefficient_value = zero
        end if
    end function ls_jj_shell_recoupling_factor

    ! Subroutine defines a unique labels for energy levels
    !
    ! Written by  G. Gaigalas                        Vilnius, May 2017 *
    ! Written by  G. Gaigalas                        Vilnius, Feb 2024 *
    subroutine assign_unique_lsj_labels(lsj_label_unit, unique_summary_unit, unique_label_unit)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: iostat_end
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use input_error_handling, only: fail_input_io
        use iso_fortran_env, only: int32
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: lsj_label_unit, unique_summary_unit, unique_label_unit

        ! Local declarations
        save  ! Preserve the legacy serial controller workspace.

        ! Local constants
        !     The maximum number of levels in the list
        integer, parameter :: lev_no  = 10000
        ! The maximum number of mixing coefficients of ASF expension
        integer, parameter :: vec_no  = 10000

        ! Local declarations
        character(len=256), pointer :: coupling(:, :)
        real(real64), pointer       :: comp(:, :), mix(:, :)
        character(len=3), pointer     :: str_no(:)
        character(len=5), pointer     :: str_j(:)
        character(len=1), pointer     :: str_p(:)
        character(len=256), pointer   :: opt_coupling(:)
        character(len=256), pointer   :: tmp_coupling(:)
        real(real64), pointer         :: energy(:), pro(:)
        real(real64), pointer         :: tmp_comp(:), tmp_mix(:)
        integer, pointer           :: mas_max(:), icount(:), iprgg(:)
        character    :: record*53
        character    :: str*7
        character(len=512) :: input_record
        real(real64) :: phase_start
        real(real64) :: max_comp
        integer      :: io_status, inum, num_lev, i1, i2, i3, i4
        integer      :: ipr, iprf, num_opt, lev_opt, new_j
        allocate(coupling(vec_no,lev_no))
        allocate(comp(vec_no,lev_no))
        allocate(mix(vec_no,lev_no))
        allocate(str_no(lev_no))
        allocate(str_j(lev_no))
        allocate(str_p(lev_no))
        allocate(tmp_coupling(vec_no))
        allocate(opt_coupling(lev_no))
        allocate(energy(lev_no))
        allocate(pro(lev_no))
        allocate(tmp_comp(vec_no))
        allocate(tmp_mix(vec_no))
        allocate(mas_max(lev_no))
        allocate(icount(lev_no))
        allocate(iprgg(lev_no))
        rewind(lsj_label_unit,iostat=io_status)
        if (io_status /= 0) call fail_input_io(io_status,lsj_label_unit,"assign_unique_lsj_labels: rewind labels")
        read (lsj_label_unit, "(1A53)", iostat=io_status) record
        if (io_status /= 0) call fail_input_io(io_status,lsj_label_unit,"assign_unique_lsj_labels: label header")
        if (len_trim(record) == 0) call fail_input_io(0,lsj_label_unit,"assign_unique_lsj_labels: empty label header")
        write(unique_label_unit, "(1A53)") record
        new_j = 0
        do while (new_j < 1)
            phase_start=profile_time()
            num_lev = 0
            new_j = 1
            do
                num_lev = num_lev + 1
                if(num_lev > lev_no) then
                    print*, "Please extand the arrays. Now we have Lev_No=",   &
                        lev_no
                    call fail_input_io(0,lsj_label_unit,"assign_unique_lsj_labels: level capacity exceeded")
                end if
                read(lsj_label_unit,"(A)",iostat=io_status) input_record
                if (io_status /= 0) call fail_input_io(io_status,lsj_label_unit,"assign_unique_lsj_labels: level record")
                if (len_trim(input_record) < 50) call fail_input_io(0,lsj_label_unit,"assign_unique_lsj_labels: short level record")
                read(input_record,"(A3,A5,5X,A1,8X,F16.9,5X,F7.3)",iostat=io_status)str_no(num_lev), &
                    str_j(num_lev),str_p(num_lev),energy(num_lev),pro(num_lev)
                if (io_status /= 0) call fail_input_io(io_status,lsj_label_unit,"assign_unique_lsj_labels: level fields")
                if (.not. ieee_is_finite(energy(num_lev)) .or. .not. ieee_is_finite(pro(num_lev))) &
                    call fail_input_io(0,lsj_label_unit,"assign_unique_lsj_labels: nonfinite level fields")
                inum = 0
                do
                    inum = inum + 1
                    if(inum > vec_no) then
                        print*, "Please extand the arrays. Now we have Vec_No=", &
                            vec_no
                        call fail_input_io(0,lsj_label_unit,"assign_unique_lsj_labels: component capacity exceeded")
                    end if
                    read(lsj_label_unit,"(A)",iostat=io_status) input_record
                    if (io_status == iostat_end) then
                        mas_max(num_lev) = inum - 1
                        exit
                    end if
                    if (io_status /= 0) call fail_input_io(io_status,lsj_label_unit,"assign_unique_lsj_labels: component record")
                    str = input_record(1:7)
                    if(str /= "       ") then
                        backspace(lsj_label_unit,iostat=io_status)
                        if (io_status /= 0) call fail_input_io(io_status,lsj_label_unit, &
                            "assign_unique_lsj_labels: next level backspace")
                        mas_max(num_lev) = inum - 1
                        exit
                    else if (len_trim(input_record) == 0) then
                        mas_max(num_lev) = inum - 1
                        new_j = 0
                        exit
                    end if
                    if (len_trim(input_record) < 37) call fail_input_io(0,lsj_label_unit, &
                        "assign_unique_lsj_labels: short component record")
                    read(input_record,"(A7,F12.8,3X,F11.8,3X,A)",iostat=io_status)str, &
                        mix(inum,num_lev),comp(inum,num_lev),coupling(inum,num_lev)
                    if (io_status /= 0) call fail_input_io(io_status,lsj_label_unit,"assign_unique_lsj_labels: component fields")
                    if (.not. ieee_is_finite(mix(inum,num_lev)) .or. .not. ieee_is_finite(comp(inum,num_lev))) &
                        call fail_input_io(0,lsj_label_unit,"assign_unique_lsj_labels: nonfinite component")
                    if (comp(inum,num_lev) < 0.0_real64) call fail_input_io(0,lsj_label_unit, &
                        "assign_unique_lsj_labels: negative component weight")
                end do
                if (io_status == iostat_end .or. new_j == 0) exit
            end do
            do i1 = 1, num_lev
                if (mas_max(i1) < 1) call fail_input_io(0,lsj_label_unit,"assign_unique_lsj_labels: level without components")
                if (comp(1,i1) <= 0.0_real64) call fail_input_io(0,lsj_label_unit, &
                    "assign_unique_lsj_labels: nonpositive leading weight")
            end do
            call profile_emit("unique_read",phase_start,items=num_lev)
            phase_start=profile_time()
            num_opt = 0
            icount = 1
            iprgg = 1
            write(unique_summary_unit,"(A)")                                                &
                "          Composition  Serial No.         Coupling"
            write(unique_summary_unit,"(A)") "                       of compos."
            write(unique_summary_unit,"(A5,A)") " J = ",trim(str_j(num_lev))
            write(unique_summary_unit,"(A)")                                                &
                "--------------------------------------------------"
            do while (num_opt < num_lev)
                max_comp = 0.0
                lev_opt = 0
                do i2 = 1, num_lev
                    if(icount(i2) == 0) cycle
                    if(comp(1,i2) > max_comp) then
                        lev_opt = i2
                        max_comp = comp(1,i2)
                    end if
                end do
                if (lev_opt == 0) call fail_input_io(0,lsj_label_unit,"assign_unique_lsj_labels: no positive remaining component")
                num_opt = num_opt + 1
                icount(lev_opt) = 0
                ipr = 1
                if(num_opt == 1) then
                    ipr = 1
                else
                    i3 = 1
                    do while (i3 < num_opt)
                        if(trim(opt_coupling(i3))==trim(coupling(ipr,lev_opt)))then
                            ipr = ipr + 1
                            if(ipr > mas_max(lev_opt)) then
                                call fail_input_io(0,lsj_label_unit,"assign_unique_lsj_labels: no unused component label")
                            end if
                            i3 = 1
                        else
                            i3 = i3 + 1
                        end if
                    end do
                end if
                iprgg(lev_opt) = ipr + iprgg(lev_opt) - 1
                if(iprgg(lev_opt) > mas_max(lev_opt)                         &
                    .and. iprgg(lev_opt) /= 1) then
                    print*,                                                   &
                        "The program is not able perform the identification for", &
                        " level = ",lev_opt
                    call fail_input_io(0,lsj_label_unit,"assign_unique_lsj_labels: no unique identification")
                end if
                if(ipr == 1) then
                    opt_coupling(num_opt) = coupling(ipr,lev_opt)
                    iprf = iprgg(lev_opt)
                    write(unique_summary_unit,"(A,I4,2X,F12.9,I5,3X,A,A,I4)")"Pos",lev_opt,     &
                        comp(ipr,lev_opt),iprf,trim(opt_coupling(num_opt))
                    if(iprf > 1) then
                        tmp_mix = 0
                        tmp_comp = 0
                        tmp_coupling = ""
                        i4 = mas_max(lev_opt)
                        tmp_mix(2:iprf)          = mix(i4-iprf+2:i4,lev_opt)
                        tmp_mix(iprf+1:i4)       = mix(2:i4-iprf+1,lev_opt)
                        mix(2:i4,lev_opt)        = tmp_mix(2:i4)
                        tmp_comp(2:iprf)         = comp(i4-iprf+2:i4,lev_opt)
                        tmp_comp(iprf+1:i4)      = comp(2:i4-iprf+1,lev_opt)
                        comp(2:i4,lev_opt)       = tmp_comp(2:i4)
                        tmp_coupling(2:iprf)    = coupling(i4-iprf+2:i4,lev_opt)
                        tmp_coupling(iprf+1:i4) = coupling(2:i4-iprf+1,lev_opt)
                        coupling(2:i4,lev_opt)  = tmp_coupling(2:i4)
                    end if
                end if
                if(ipr > 1) then
                    tmp_mix = 0
                    tmp_comp = 0
                    tmp_coupling = ""
                    i4 = mas_max(lev_opt)
                    tmp_mix(1:1) = mix(ipr:ipr,lev_opt)
                    tmp_mix(i4-ipr+2:i4)      = mix(1:ipr-1,lev_opt)
                    tmp_mix(2:i4-ipr+1)       = mix(ipr+1:i4,lev_opt)
                    mix(1:i4,lev_opt)         = tmp_mix(1:i4)
                    tmp_comp(1:1)             = comp(ipr:ipr,lev_opt)
                    tmp_comp(i4-ipr+2:i4)     = comp(1:ipr-1,lev_opt)
                    tmp_comp(2:i4-ipr+1)      = comp(ipr+1:i4,lev_opt)
                    comp(1:i4,lev_opt)        = tmp_comp(1:i4)
                    tmp_coupling(1:1)         = coupling(ipr:ipr,lev_opt)
                    tmp_coupling(i4-ipr+2:i4) = coupling(1:ipr-1,lev_opt)
                    tmp_coupling(2:i4-ipr+1)  = coupling(ipr+1:i4,lev_opt)
                    coupling(1:i4,lev_opt)    = tmp_coupling(1:i4)
                    icount(lev_opt) = 1
                    num_opt = num_opt - 1
                end if
            end do
            write(unique_summary_unit,"(A)")                                                &
                "--------------------------------------------------"
            write(unique_summary_unit,"(A)")""
            call profile_emit("unique_assign_summary",phase_start,items=num_lev)
            phase_start=profile_time()
            do i1 = 1, num_lev
                write(unique_label_unit,"(A3,A5,5X,A1,8X,F16.9,5X,F7.3,A)")                 &
                    str_no(i1),str_j(i1),str_p(i1),energy(i1),pro(i1),"%"
                do i2 = 1, mas_max(i1)
                    write(unique_label_unit,"(A7,F12.8,3X,F11.8,3X,A)")                       &
                        str,mix(i2,i1),comp(i2,i1),trim(coupling(i2,i1))
                end do
            end do
            if(new_j < 1) write(unique_label_unit,"(A)") " "
            call profile_emit("unique_write",phase_start,items=num_lev)
        end do
        deallocate(coupling,comp,mix,str_no,str_j,str_p,tmp_coupling,opt_coupling, &
            energy,pro,tmp_comp,tmp_mix,mas_max,icount,iprgg)
    end subroutine assign_unique_lsj_labels
end module jj2lsj_conversion
