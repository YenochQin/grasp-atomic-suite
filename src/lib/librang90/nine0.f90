!*******************************************************************
!                                                                  *
      subroutine nine0(j1,j2,j3,l1,l2,l3,k1,k2,k3,aa)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
!                                                                  *
!     |  J1/2  J2/2  J3/2 |                                        *
!     |  L1/2  L2/2  L3/2 |                                        *
!     |  K1/2  K2/2    0  |                                        *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                             March    1995  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use cons_c,          only: zero, one
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use sixj_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j1, j2, j3, l1, l2, l3, k1, k2, k3
      real(double), intent(out) :: aa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: ifa
      real(double) :: a, b
!-----------------------------------------------
      if (j1 == 0) then
         call sixj (l2, k2, j2, k3, l3, l1, 0, a)
         b = dble((j2 + 1)*(l1 + 1))
         ifa = k2 + j2 + l3 + l1
      else if (j2 == 0) then
         call sixj (l3, k3, j3, k1, l1, l2, 0, a)
         b = dble((j3 + 1)*(l2 + 1))
         ifa = k3 + j3 + l1 + l2
      else if (j3 == 0) then
         call sixj (l1, k1, j1, k2, l2, l3, 0, a)
         b = dble((j1 + 1)*(l3 + 1))
         ifa = k1 + j1 + l2 + l3
      else if (l1 == 0) then
         call sixj (k2, j2, l2, j3, k3, k1, 0, a)
         b = dble((l2 + 1)*(k1 + 1))
         ifa = j2 + l2 + k3 + k1
      else if (l2 == 0) then
         call sixj (k3, j3, l3, j1, k1, k2, 0, a)
         b = dble((l3 + 1)*(k2 + 1))
         ifa = j3 + l3 + k1 + k2
      else if (l3 == 0) then
         call sixj (k1, j1, l1, j2, k2, k3, 0, a)
         b = dble((l1 + 1)*(k3 + 1))
         ifa = j1 + l1 + k2 + k3
      else if (k1 == 0) then
         call sixj (j2, j3, j1, l3, l2, k2, 0, a)
         b = dble((j1 + 1)*(k2 + 1))
         ifa = j3 + j1 + l2 + k2
      else if (k2 == 0) then
         call sixj (j3, j1, j2, l1, l3, k3, 0, a)
         b = dble((j2 + 1)*(k3 + 1))
         ifa = j1 + j2 + l3 + k3
      else if (k3 == 0) then
         call sixj (j1, j2, j3, l2, l1, k1, 0, a)
         b = dble((j3 + 1)*(k1 + 1))
         ifa = j2 + j3 + l1 + k1
      else
         a = zero
         b = one
      endif
      aa = a/dsqrt(b)
      if (mod(ifa,4) /= 0) aa = -aa
      return
      end subroutine nine0
