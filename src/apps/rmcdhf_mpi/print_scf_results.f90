module print_scf_results_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: print_scf_results
contains

    !***********************************************************************
    !                                                                      *
    subroutine print_scf_results
        ! Module dependencies
        !                                                                      *
        !   This subroutine is now called only in GETOLD to print out all
        !   orbitals so that the user can do cut-and-paste in supplying the
        !   orbitals to be varied. I.e., LFIX(LOC) will always be FALSE and
        !   are thus commented out.
        ! Xinghong HE 6 Apr 1998
        !                                                                      *
        !   Call(s) to: [LIB92]: CONVRT.                                       *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        use orbital_csf_basis, only: orbital_count, orbital_kappa, orbital_principal_n, orbital_symmetry_labels
        use standard_io_units, only: stderr_unit
        use integer_format, only: format_integer
        implicit none (type, external)

        ! Local declarations
        integer :: iend, orbital_index, loc, ibeg, text_length
        character :: record*80, number_text*2
        !-----------------------------------------------
        !
        !
        !   Print the list of all subshell radial wavefunctions
        !
        iend = 0
        do orbital_index = 1, orbital_count
            !         LOC = IORDER(I)
            loc = orbital_index
            ! For the commenting-out of the IF...ENDIF, see comments in the header
            !         IF (.NOT. LFIX(LOC)) THEN
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
            call format_integer (orbital_principal_n(loc), number_text, text_length)
            iend = ibeg + text_length - 1
            record(ibeg:iend) = number_text(1:text_length)
            if (orbital_kappa(loc) < 0) then
                text_length = 1
            else
                text_length = 2
            end if
            ibeg = iend + 1
            iend = ibeg + text_length - 1
            record(ibeg:iend) = orbital_symmetry_labels(loc)(1:text_length)
            !         ENDIF
        end do
        if (iend > 1) write (stderr_unit, *) record(1:iend)
        !
        return
    end subroutine print_scf_results
end module print_scf_results_interface
