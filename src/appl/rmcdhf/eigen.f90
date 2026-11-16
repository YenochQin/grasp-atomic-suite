!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function eigen (j)
!                                                                      *
!   This function computes an estimate of the energy of orbital J .    *
!                                                                      *
!   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use def_c
      use grid_c
      use orb_c
      use pote_c
      use tatb_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use quad_i
      use dpbdt_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(double) :: piece1, piece2, piece3, piece4, piece5
!-----------------------------------------------
!
!
!   Initialization
!
      mtp = mf(j)
!
!   Exchange term
!
      ta(1) = 0.0d00
      do i = 2, mtp
         ta(i) = (pf(i,j)*xq(i)-qf(i,j)*xp(i))*rpor(i)
      end do
      call quad (piece1)
      piece1 = c*piece1
!
!   Direct term
!
      ta(1) = 0.0d00
      do i = 2, mtp
         ta(i) = (pf(i,j)**2 + qf(i,j)**2)*yp(i)*rpor(i)
      end do
      call quad (piece2)
!
!   Kinetic energy terms
!
      ta(1) = 0.0d00
      do i = 2, mtp
         ta(i) = qf(i,j)**2*rp(i)
      end do
      call quad (piece3)
      piece3 = 2.0d00*c*c*piece3
!
      ta(1) = 0.0d00
      do i = 2, mtp
         ta(i) = (pf(i,j)*qf(i,j))*rpor(i)
      end do
      call quad (piece4)
      piece4 = -2.0d00*dble(nak(j))*c*piece4
!
      call dpbdt (j)
      ta(1) = 0.0d00
      do i = 2, mtp
         ta(i) = pf(i,j)*tb(i) - qf(i,j)*ta(i)
      end do
      call quad (piece5)
      piece5 = c*piece5/h
!
!   Assembly
!
      eigen = piece1 + piece2 + piece3 + piece4 + piece5
!
      return
!
      end function eigen
