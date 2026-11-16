!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function rinti (j, k, mode)
!                                                                      *
!   The value of this  function is the one-electron integral I (J,K)   *
!   for  orbitals  J, K. The analytical expression for this quantity   *
!   is given as  eq (9) in  I P Grant, B J McKenzie, P H Norrington,   *
!   D F Mayers, and N C Pyper,  Computer  Phys Commun 21 (1980) 211 .  *
!                                                                      *
!   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:23   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use debug_c
      use def_c,           only: c
      use grid_c
      use npot_c,          only: zz
      use orb_c
      use tatb_c,          only: ta, tb, mtp
      use wave_c,          only: mf,pf,qf
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dpbdt_i
      use quad_i
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
      integer  :: k
      integer, intent(in) :: mode
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(double) :: piece1, piece2, piece3, piece4
!-----------------------------------------------
!
!   Stop if orbitals J and K have different kappa values
!
      if (nak(j) /= nak(k)) then
         write (*, 300) np(j), nh(j), np(k), nh(k)
         stop
      endif
!
      mtp = max(mf(j),mf(k))
!
!   Kinetic energy contribution
!
!   Piece involving derivatives
!
      call dpbdt (k)
      ta(1) = 0.d0
      do i = 2, mtp
         ta(i) = qf(i,j)*ta(i) - pf(i,j)*tb(i)
      end do
      call quad (piece1)
      piece1 = c*piece1/h
!
!   Pieces not involving derivatives
!
      ta(1) = 0.d0
      do i = 2, mtp
         ta(i) = rp(i)*qf(i,j)*qf(i,k)
      end do
      call quad (piece2)
      piece2 = -2.d0*c*c*piece2
!
      ta(1) = 0.d0
      do i = 2, mtp
         ta(i) = rpor(i)*(pf(i,j)*qf(i,k) + qf(i,j)*pf(i,k))
      end do
      call quad (piece3)
      piece3 = piece3*c*dble(nak(k))
!
!   Contribution from nuclear potential only if MODE is 0
!
      if (mode == 0) then
!
         ta(1) = 0.d0
         do i = 2, mtp
            ta(i) = rpor(i)*zz(i)*(pf(i,j)*pf(i,k) + qf(i,j)*qf(i,k))
         end do
         call quad (piece4)
         piece4 = -piece4
!
      else
         piece4 = 0.d0
      endif
!
      rinti = piece1 + piece2 + piece3 + piece4
!
!   Debug printout
!
      if (mode==0 .and. ldbpr(4)) write (99, 301) np(j), nh(j), np(k), nh(k), &
         rinti
      if (mode/=0 .and. ldbpr(5)) write (99, 302) np(j), nh(j), np(k), nh(k), &
         rinti
!
      return
!
  300 format('RINTI: Attempt to calculate I(',1i2,1a2,',',1i2,1a2,')')
  301 format(/,' I (',1i2,1a2,',',1i2,1a2,') = ',1p,d19.12,/)
  302 format(/,' K (',1i2,1a2,',',1i2,1a2,') = ',1p,d19.12,/)
      return
!
      end function rinti
