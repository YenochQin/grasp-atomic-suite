!***********************************************************************
!                                                                      *
      subroutine qed(jstate, slfint, ucf)
!                                                                      *
!   This  routine estimates corrections to the  energy levels due to   *
!   self-energy.                                                       *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, IQ, QUAD, RALLOC, SCREEN.      *
!               [RCI92]: FZALF, HOVLAP.                                *
!                                                                      *
!                                           Last update: 30 Oct 1992   *
!   Modified by Xinghong He                 Last update: 24 Jun 1997   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: nnnw, nnnp
      use def_c
      use eigv_c
      use npar_c
      use orb_c, only: nw, ncf, np, nak
      use prnt_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use iq_i
      use ratden_i
      use fzalf_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jstate
      real(double), intent(out) :: slfint(nnnw)
      real(double), intent(out) :: ucf(1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: maxiter, j, i, ii, npj, kappa, mfj
      real(double), dimension(nnnp) :: ptemp, qtemp
      real(double) :: ucfj, zeff, ratio, valu
      character :: npchar, nakchar*2
!-----------------------------------------------
!
! Pre-set tolerable number for iteration in finding effective
! nuclear charge.
!
      maxiter = 20
!
! Modified so that UCFJ describes the current eigenstate
!
      do j = 1, nw
         ucfj = 0.0d00
!         DO 3 I = 1,NVEC
         i = jstate
         do ii = 1, ncf
            ucfj = ucfj + dble(iq(j,ii))*evec(ii + (i - 1)*ncf)**2
         end do
!    3    CONTINUE
!     print *, ucfj,'ucf'
         ucf(j) = ucfj
! zou    UCF(J) = UCFJ/DBLE (NCF)
      end do
!
      do j = 1, nw
!
         npj = np(j)
!
         if (npj <= 8) then
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
!
            zeff = z
            ratio = ratden(ptemp,qtemp,mfj,npj,kappa,zeff)
            valu = ratio*fzalf(npj,kappa,zeff)/dble(npj**3)
!
            slfint(j) = valu*zeff**4/(pi*c**3)
!        print *, 'No. orb.=',j,' Zeff = ',zeff
!    &   , 'Scale= ',ratio
!    &   , 'S.E. = ',slfint(j)*2*13.6058,slfint(j)/ratio*2*13.6058
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
      end subroutine qed


!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function ratden (p, q, mtpo, np, kappa, z)
!                                                                      *
!   This subprogram computes the overlap of the orbital tabulated in   *
!   the arrays  P  and  Q  with maximum tabulation point  MTPO  with   *
!   a hydrogenic orbital with parameters  NP  KAPPA  Z .               *
!                                                                      *
!   Call(s) to: [LIB92]: DCBSRW, QUAD.                                 *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Switches:
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: nnnp
      use grid_c
      use horb_c, only: ph, qh
      use tatb_c, only: mtp, ta
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dcbsrw_i
      use quad_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: mtpo
      integer  :: np
      integer  :: kappa
      real(double)  :: z
      real(double), intent(in) :: p(nnnp)
      real(double), intent(in) :: q(nnnp)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mtph, i, k
      real(double) :: eh, pzh, result, result1
!-----------------------------------------------
!
!
!
!   Set up the hydrogenic orbital
!
      call dcbsrw (np, kappa, z, eh, pzh, ph, qh, mtph)
!
!   Compute the overlap
!
      mtp = min(mtph,mtpo)
      do i = 2, mtp
         if (rp(i) > 0.0219) cycle
         k = i
      end do
      mtp = k
      ta(1) = 0.0d00
      ta(2:mtp) = (p(2:mtp)*p(2:mtp)+q(2:mtp)*q(2:mtp))*rp(2:mtp)
!         TA(I) = (P(I)*PH(I)+Q(I)*QH(I))*RP(I)
      call quad (result)
      ta(2:mtp) = (ph(2:mtp)*ph(2:mtp)+qh(2:mtp)*qh(2:mtp))*rp(2:mtp)
      call quad (result1)
!
      ratden = result/result1
!
      return
      end function ratden
