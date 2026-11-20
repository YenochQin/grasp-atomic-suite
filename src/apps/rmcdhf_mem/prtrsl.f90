!***********************************************************************
!                                                                      *
      subroutine prtrsl
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_kappa, orbital_principal_n, &
          orbital_symmetry_labels
      use fixed_orbital_flags, only: fixed_orbital_count, is_orbital_fixed, max_orbitals
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use integer_format, only: convrt
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
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
         call convrt (orbital_principal_n(loc), number_text, text_length)
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
      end subroutine prtrsl
