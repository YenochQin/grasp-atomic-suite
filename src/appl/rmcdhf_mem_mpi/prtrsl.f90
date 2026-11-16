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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use orb_c
      use fixd_c
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iend, i, loc, ibeg, lenth
      character :: record*80, cnum*2
!-----------------------------------------------
!
!
!   Print the list of all subshell radial wavefunctions
!
      iend = 0
      do i = 1, orbital_count
!         LOC = IORDER(I)
         loc = i
! For the commenting-out of the IF...ENDIF, see comments in the header
!         IF (.NOT. LFIX(LOC)) THEN
         if (iend > 75) then
            write (stderr_unit, *) record(1:iend)
            iend = 0
         endif
         if (iend > 0) then
            ibeg = iend + 1
            iend = ibeg
            record(ibeg:iend) = ' '
         endif
         ibeg = iend + 1
         call convrt (orbital_principal_n(loc), cnum, lenth)
         iend = ibeg + lenth - 1
         record(ibeg:iend) = cnum(1:lenth)
         if (orbital_kappa(loc) < 0) then
            lenth = 1
         else
            lenth = 2
         endif
         ibeg = iend + 1
         iend = ibeg + lenth - 1
         record(ibeg:iend) = orbital_symmetry_labels(loc)(1:lenth)
!         ENDIF
      end do
      if (iend > 1) write (stderr_unit, *) record(1:iend)
!
      return
      end subroutine prtrsl
