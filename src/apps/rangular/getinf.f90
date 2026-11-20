!***********************************************************************
!                                                                      *
      subroutine getinf
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use numeric_constants, only: eleven, eps, four, half, one, seven, ten, tenth, three, two, zero
      use csf_block_cutoffs, only: iccut, iccutblk
      use calculation_input_defaults, only: inpci, ndef, ndump
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yes_no_input, only: getyn
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
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
            iccut(i) = 0
         end do
      else
         if (ndef /= 0) then
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
               read *, iccut(i)
               write(angular_log_unit,*) iccut(i), "! ICCUT FOR BLOCK",i
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
      end subroutine getinf
