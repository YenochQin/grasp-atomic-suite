!***********************************************************************
!                                                                      *
      subroutine qed_slfen(slfint)
!                                                                      *
!   This  routine estimates the F(Z\alpha) function of self energy for *
!   each orbital.                                                      *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, IQ, QUAD, RALLOC, SHIELD.      *
!               [RCI92]: FZALF, HOVLAP.                                *
!                                                                      *
!   Modified from subroutine QED by Yu Zou, Last update: 13 Mar 2000   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use parameter_def,   only: nnnw, nnnp
      use def_c
      use eigv_c
      use npar_c
      use orb_c, only: nw, np, nak
      use prnt_c
      use wave_c
      use qedcut_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ratden_i
      use fzalf_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double), dimension(nnnw), intent(out) :: slfint
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: maxiter, j, npj, kappa, mfj, i, npjmax
!GG      REAL(DOUBLE), DIMENSION(1) :: UCF
      real(double), dimension(nnnp) :: ptemp, qtemp
      real(double) :: zeff, ratio, valu
      character :: npchar, nakchar*2
!-----------------------------------------------
!
! Pre-set tolerable number for iteration in finding effective
! nuclear charge.
!
      maxiter = 20
!Per
      if (nqedcut.eq.1) then
         npjmax = nqedmax
      else
         npjmax = 8
      end if
!Per
!
      do j = 1, nw
!
         npj = np(j)
!
         if (npj <= npjmax) then
!
!   Only orbitals with principal quantum number 8 or less can
!   be treated by this section of code
!
            kappa = nak(j)
!
!   Begin by transferring the function to a temporary array
!
            mfj = mf(j)
!
            ptemp(1) = 0.0d00
            qtemp(1) = 0.0d00
            do i = 2, mfj
               ptemp(i) = pf(i,j)
               qtemp(i) = qf(i,j)
            end do
            zeff = z
            ratio = ratden(ptemp,qtemp,mfj,npj,kappa,zeff)
            valu = ratio*fzalf(npj,kappa,zeff)/dble(npj**3)
            slfint(j) = valu*zeff**4/(pi*c**3)
!
         else
!
!   The self-energy for orbitals with principal quantum number
!   greater than 8 is set to zero
!
            slfint(j) = 0.0d00
!
         endif
!
      end do
!
!   Deallocate storage for the `generalised occupation numbers'
!
!
      return
      end subroutine qed_slfen
