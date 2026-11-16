!*******************************************************************
!                                                                  *
      subroutine nine(j1,j2,j3,l1,l2,l3,k1,k2,k3,i,inn,aa)
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
!                                                                  *
!     |  J1/2  J2/2  J3/2 |                                        *
!     |  L1/2  L2/2  L3/2 |                                        *
!     |  K1/2  K2/2  K3/2 |                                        *
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
      use, intrinsic :: iso_fortran_env, only: real64
      use cons_c,          only: zero
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ittk_i
      use nine0_i
      use sixj_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j1, j2, j3, l1, l2, l3, k1, k2, k3
      integer, intent(in)  :: i
      integer, intent(out) :: inn
      real(real64)  :: aa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: n1, n2, n3, n4, n5, n6, max_, min_, ix
      real(real64) :: s1, s2, s3, x
!-----------------------------------------------
      if (i == 1) then
         inn = 0
         if (ittk(j1,j2,j3) == 0) return
         if (ittk(l1,l2,l3) == 0) return
         if (ittk(k1,k2,k3) == 0) return
         if (ittk(j1,l1,k1) == 0) return
         if (ittk(j2,l2,k2) == 0) return
         if (ittk(j3,l3,k3) == 0) return
         inn = 1
         return
      endif
      if (j1*j2*j3*l1*l2*l3*k1*k2*k3 == 0) then
         inn = 1
         call nine0 (j1, j2, j3, l1, l2, l3, k1, k2, k3, aa)
      else
         n1 = iabs(j1 - k3)
         n2 = iabs(l3 - j2)
         n3 = iabs(l1 - k2)
         n4 = iabs(j2 - l3)
         n5 = iabs(k2 - l1)
         n6 = iabs(j1 - k3)
         max_ = max0(n1,n2,n3,n4,n5,n6)
         n1 = j1 + k3
         n2 = l3 + j2
         n3 = j2 + l3
         n4 = k2 + l1
         n5 = j1 + k3
         n6 = l1 + k2
         min_ = min0(n1,n2,n3,n4,n5,n6)
         inn = 1
         aa = zero
         do ix = max_, min_, 2
            call sixj (j1, j2, j3, l3, k3, ix, 0, s1)
            call sixj (l1, l2, l3, j2, ix, k2, 0, s2)
            call sixj (k1, k2, k3, ix, j1, l1, 0, s3)
            x = s1*s2*s3*dble(ix + 1)
            if (mod(ix,2) /= 0) x = -x
            aa = x + aa
         end do
      endif
      return
      end subroutine nine
