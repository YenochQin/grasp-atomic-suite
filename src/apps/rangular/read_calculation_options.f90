module read_calculation_options_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_calculation_options
contains

    !***********************************************************************
    !                                                                      *
    subroutine read_calculation_options
        ! Module dependencies
        use calculation_file_units, only: angular_log_unit
        !                                                                      *
        !   Interactively determines data governing the generation of MCP co-  *
        !   efficients.                                                        *
        !                                                                      *
        !   Call(s) to: [LIB92]: CONVRT, GETYN.                                *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 19 Dec 1992   *
        !   Modified Xinghong He                  Last revision: 30 Jun 1998   *
        !                                                                      *
        !   File shared (hard link) by mcpmpi, mcpblk
        !
        !   Updated to treat ICCUT for block
        !
        !***********************************************************************
        use csf_block_cutoffs, only: active_csf_cutoffs
        use calculation_input_defaults, only: nondefault_options_flag
        use standard_io_units, only: stderr_unit
        use mcp_generation_options, only: diagonal_only, use_first_order_restriction
        use hamiltonian_symmetry_blocks, only: symmetry_block_count
        use yes_no_input, only: read_yes_no
        implicit none (type, external)

        ! Local declarations
        integer :: i
        !-----------------------------------------------

        !-----------------------------------------------------------------------
        !
        !   Determine the physical effects specifications
        !
        !     Commenting out the EAL option
        !     IF (NDEF .NE. 0) THEN
        !        WRITE  (istde,*) 'Generate MCP coefficients only for'
        !    & , ' diagonal matrix elements? '
        !        WRITE (istde,*) '(This is appropriate to (E)AL calculation):'
        !        DIAG = GETYN ()
        !     ELSE
        !        DIAG = .FALSE.
        !     ENDIF
        diagonal_only = .false.
        if (diagonal_only) then
            use_first_order_restriction = .false.
            do i = 1,100
                active_csf_cutoffs(i) = 0
            end do
        else
            if (nondefault_options_flag /= 0) then
                !            WRITE (istde,*) 'Treat contributions of some CSFs', &
                !                            ' as first-order perturbations?'
                !            LFORDR = GETYN ()
                use_first_order_restriction = .true.
            else
                use_first_order_restriction = .false.
            end if
            if (use_first_order_restriction) then
                write (stderr_unit,*) "The contribution of CSFs 1 -- ICCUT will",&
                    " be treated variationally;"
                write (stderr_unit,*) "the remainder perturbatively; enter ICCUT:"
                do i = 1,symmetry_block_count
                    write(stderr_unit,*) "Give ICCUT for block",i
                    read *, active_csf_cutoffs(i)
                    write(angular_log_unit,*) active_csf_cutoffs(i), "! ICCUT FOR BLOCK",i
                    !    1          READ *, ICCUT(i)
                    !               IF ((ICCUT(i).LE.1).OR.(ICCUT(i).GE.ncfblk(i))) THEN
                    !                  WRITE (istde,*) 'GETINF: ICCUT must be greater than 1', &
                    !                                  ' and less than ',ncfblk(i)
                    !                  WRITE (istde,*) ' please reenter ICCUT:'
                    !                  GOTO 1
                    !               ENDIF
                end do
            end if
        end if

        return
    end subroutine read_calculation_options
end module read_calculation_options_interface
