!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function slater (ia, ib, ic, id, k)
!                                                                      *
!   The value of this  function is the Slater integral                 *
!                                                                      *
!                               k                                      *
!                              R (abcd)                                *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD, YZK.                                    *
!                                                                      *
!                                         Last revision: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:44   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use debug_c
      use grid_c
      use orb_c
      use tatb_c,          only: ta, tb, mtp
      use wave_c,          only: mf, qf, pf
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yzk_i
      use quad_i
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: ia
      integer :: ib
      integer :: ic
      integer :: id
      integer :: k
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(double) :: result
!-----------------------------------------------
!
!
      call yzk (k, ib, id)
!
!   Multiply by second term, and obtain result by integration
!
      if (k==0 .and. ib==id) then
         mtp = min(mf(ia),mf(ic))
      else
         mtp = min(min(mtp,mf(ia)),mf(ic))
      endif
!
      ta(1) = 0.0d00
      do i = 2, mtp
         ta(i) = (pf(i,ia)*pf(i,ic) + qf(i,ia)*qf(i,ic))*rpor(i)*tb(i)
      end do
!
      call quad (result)
      slater = result
!
!   Debug printout
!
      if (ldbpr(10)) write (99, 300) k, np(ia), nh(ia), np(ib), nh(ib), np(ic)&
         , nh(ic), np(id), nh(id), slater
!
      return
!
  300 format(/,'  (',1i1,')',/,' R   (',1i2,1a2,',',1i2,1a2,';',1i2,1a2,',',1i2&
         ,1a2,') ','= ',1p,d19.12)
      return
!
      end function slater
