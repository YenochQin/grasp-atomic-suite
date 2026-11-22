module print_uninitialized_orbitals_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: print_uninitialized_orbitals
contains

    !   Prints a list of subshells that remain to be estimated.  ALL  is
    !   .TRUE.  if all subshells have been estimated and  .FALSE. other-
    !   wise.
    !
    ! Written by Farid A. Parpia. Last revision: 15 Dec 1992
    subroutine print_uninitialized_orbitals(all)
        ! Module dependencies
        use initial_orbital_selection, only: orbital_initialized
        use orbital_csf_basis, only: orbital_count, orbital_principal_n, orbital_kappa, orbital_symmetry_labels
        use standard_io_units, only: stderr_unit
        use integer_format, only: format_integer
        implicit none (type, external)

        ! Arguments
        logical , intent(out) :: all

        !  Local declarations
        ! ...  /ORB2/
        !      COMMON /ORB2/ NCF, NW, PNTRIQ
        !      REAL   PNTRIQ
        integer :: orbital_index, iend, ibeg, text_length
        character(len=12) :: number_text
        character(len=80) :: record*80
        ! Determine if there are any subshell radial wavefunctions that
        ! remain to be estimated
        all = .true.
        do orbital_index = 1, orbital_count
            if (orbital_initialized(orbital_index)) cycle
            all = .false.
            exit
        end do
        if (all) return
        ! Print a list of subshell radial wavefunctions that remain to
        ! be estimated; this list is no more than 80 characters wide
        write (stderr_unit, *) "The following subshell radial wavefunctions ", &
            "remain to be estimated:"
        iend = 0
        do orbital_index = 1, orbital_count
            if (orbital_initialized(orbital_index)) cycle
            if (iend > 75) then
                write (stderr_unit, *) record(1:iend)
                iend = 0
            end if
            if (iend > 0) then
                ibeg = iend + 1
                iend = ibeg
                record(ibeg:iend) = " "
            end if
            ibeg = iend + 1
            call format_integer (orbital_principal_n(orbital_index), number_text, text_length)
            iend = ibeg + text_length - 1
            record(ibeg:iend) = number_text(1:text_length)
            if (orbital_kappa(orbital_index) < 0) then
                text_length = 1
            else
                text_length = 2
            end if
            ibeg = iend + 1
            iend = ibeg + text_length - 1
            record(ibeg:iend) = orbital_symmetry_labels(orbital_index)(1:text_length)
        end do
        if (iend > 1) write (stderr_unit, *) record(1:iend)
        return
    end subroutine print_uninitialized_orbitals
end module print_uninitialized_orbitals_interface
