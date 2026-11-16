!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function skint (rac, ia, ic, rbd, ib, id, k, iw)
!                                                                      *
!   This routine evaluates transverse interaction integrals:           *
!                                                                      *
!                          (k)                                         *
!                         S   (a,c;b,d;w)                              *
!                                                                      *
!   where w = wac if IW = 1, and w= wbd if IW = 2.                     *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!               [RCI92]: ZKF.                                          *
!                                           Last update: 06 Nov 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: nnnp
      use bess_c
      use debug_c
      use grid_c, only: r, rp, rpor
      use orb_c
      use tatb_c, only: mtp, ta, tb
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use zkf_i
      use quad_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ia, ib, ic, id
      integer, intent(in) :: k, iw
      real(double), dimension(nnnp), intent(in) :: rac
      real(double), dimension(nnnp), intent(in) :: rbd
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mxrac, mxrbd, i
      real(double), dimension(nnnp) :: tkeep
      real(double) :: epsi, w, wk, valu
!-----------------------------------------------
!
!
      epsi = 1.0d-10
!
      w = wij(iw)
      wk = dble(k + k + 1)
      mxrac = min(mf(ia),mf(ic))
      mxrbd = min(mf(ib),mf(id))
!
!           (k-1)
!  Compute Z     (rho  ; s)
!                    ac
!
      ta(:mxrac) = rac(:mxrac)
      mtp = mxrac
      call zkf (k - 1, ia, ic)
!
      if (dabs(w) < epsi) then
!
!   W = 0 case
!
         tkeep(:mtp) = tb(:mtp)
         call zkf (k + 1, ia, ic)
         mtp = min(mtp,mxrbd)
         ta(1) = 0.0d00
         ta(2:mtp) = rbd(2:mtp)*rpor(2:mtp)*(tb(2:mtp)-tkeep(2:mtp))
         call quad (valu)
         skint = wk*valu*0.5d00
!
      else
!
!   Finite w: see I P Grant and B J McKenzie, J Phys B: At Mol Phys,
!   13 (1980) 2671-2681
!
         tkeep(:mtp) = tb(:mtp)
         ta(:mtp) = -ta(:mtp)*bessj(1,iw,:mtp)
         call zkf (k - 1, ia, ic)
!
         mtp = min(mtp,mxrbd)
         ta(1) = 0.0d00
         ta(2:mtp) = ((1.0d00 + bessn(2,iw,2:mtp))*tb(2:mtp)-tkeep(2:mtp)*bessn&
            (2,iw,2:mtp))*rbd(2:mtp)/r(2:mtp)**2*rpor(2:mtp)
         call quad (valu)
         skint = (wk/w)**2*valu
         ta(:mxrbd) = rbd(:mxrbd)*(1.0d00 + bessj(2,iw,:mxrbd))
         mtp = mxrbd
         call zkf (k + 1, ib, id)
         mtp = min(mtp,mxrac)
         ta(1) = 0.0d00
         ta(2:mtp) = rac(2:mtp)*(1.0d00 + bessn(1,iw,2:mtp))*tb(2:mtp)*r(2:mtp)&
            *rp(2:mtp)
         call quad (valu)
         skint = skint - valu*w*w/dble((2*k + 3)*(2*k - 1))
!
      endif
!
      if (ldbpr(11)) write (99, 300) k, np(ia), nh(ia), np(ic), nh(ic), np(ib)&
         , nh(ib), np(id), nh(id), iw, valu
!
      return
!
  300 format('  (',1i2,')'/,'S     (',1i2,1a2,',',1i2,1a2,'|',1i2,1a2,',',1i2,&
         1a2,';',1i2,') = ',1p,d19.12)
      return
!
      end function skint
