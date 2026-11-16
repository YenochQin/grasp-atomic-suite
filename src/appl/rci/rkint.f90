!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function rkint (rac, ia, ic, rbd, ib, id, k, iw)
!                                                                      *
!   This routine evaluates the transverse interaction integrals.  If   *
!   IW = 0, it calulates the U(r1,r2) integral; otherwise, it calcu-   *
!   lates  R bar (k; a c | b d ; w)  with  w = wac  if  IW = 1,  and   *
!   w = wbd  if IW = 2.                                                *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!               [RCI92]: ZKF.                                          *
!                                           Last update: 15 Oct 1992   *
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
      use grid_c, only: rpor
      use orb_c
      use tatb_c,  only: mtp, ta, tb
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
      integer  :: ia, ic
      integer, intent(in) :: ib, id
      integer  :: k
      integer, intent(in) :: iw
      real(double), dimension(nnnp), intent(in) :: rac
      real(double), dimension(nnnp), intent(in) :: rbd
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mxrbd, mxrac, i
      real(double) :: result
!-----------------------------------------------
!
      mxrbd = min (mf(ib),mf(id))
      mxrac = min (mf(ia),mf(ic))
!
      if (iw == 0) then
!
!   IW = 0
!
         ta(:mxrac) = rac(:mxrac)
         mtp = mxrac
         call zkf (k, ia, ic)
         mtp = min(mtp,mxrbd)
         ta(1) = 0.0d00
         ta(2:mtp) = rbd(2:mtp)*tb(2:mtp)*rpor(2:mtp)
!
      else
!
!   IW = 1,2
!
         ta(:mxrac) = rac(:mxrac)*(1.0d00 + bessj(1,iw,:mxrac))
         mtp = mxrac
         call zkf (k, ia, ic)
         mtp = min(mtp,mxrbd)
         ta(1) = 0.0d00
         ta(2:mtp) = rbd(2:mtp)*(1.0d00 + bessn(1,iw,2:mtp))*tb(2:mtp)*rpor(2:&
            mtp)
!
      endif
!
      call quad (result)
      rkint = result
!
!   Debug printout if option set
!
      if (ldbpr(11)) write (99, 300) k, np(ia), nh(ia), np(ic), nh(ic), np(ib)&
         , nh(ib), np(id), nh(id), iw, result
!
      return
!
  300 format('_ (',1i2,')'/,'R     (',1i2,1a2,',',1i2,1a2,'|',1i2,1a2,',',1i2,&
         1a2,';',1i2,') = ',1p,d19.12)
      return
!
      end function rkint
