module sixj_i
    implicit none (type, external)
    private
    public :: sixj
contains
!*******************************************************************
!                                                                  *
      subroutine sixj(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6,check_triangle_conditions,six_j_value)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | I/2  J/2  K/2 |                                            *
!     | L/2  M/2  N/2 |          (29.1A) [J.B.77]                  *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vanderbilt University,  Nashville               October  1996  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use cons_c,          only: zero, one
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ixjtik_i
      use sixj5_i
      use sixj1_i
      use sixj35_i
      use sixj2_i
!      USE gracah1_I
      use dracah_i
      use sixj3_i
      use sixj4_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: two_j1, two_j2, two_j3, two_j4, two_j5, two_j6
      integer, intent(in) :: check_triangle_conditions
      real(real64), intent(out) :: six_j_value
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ifa
      real(real64) :: small_six_j_cache(0:4, 0:4, 0:4, 0:4, 0:4, 0:4)
      real(real64) :: uncached_six_j_sentinel, a
      logical :: store_in_six_j_cache
!     Persistent memoization is intentional (DATA already implies SAVE).
!     State remains per-thread when this source is compiled with OpenMP.
      save :: small_six_j_cache
!-----------------------------------------------
      data small_six_j_cache/ 15625*1.e-20_real64/
      data uncached_six_j_sentinel/ 1.e-20_real64/
!     The jj2lsj OpenMP object build gives each thread its own lazy cache.
!     Ordinary suite libraries compile this sentinel as a comment.
!$OMP THREADPRIVATE(small_six_j_cache)
      six_j_value = zero
      if (check_triangle_conditions /= 0) then
!
!     CHESKED TRIANGULAR CONDITIONS
!
         if (ixjtik(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6) == 0) return
      endif
      store_in_six_j_cache = .false.
      if (max0(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6) <= 4) then
         six_j_value = small_six_j_cache(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6)
         if (six_j_value == uncached_six_j_sentinel) then
            store_in_six_j_cache = .true.
         else
            return
         endif
      endif
!
!     CALCULATED IN CASE WHEN ONE OF PERAMETERS EQUAL 0
!
      if (two_j1*two_j2*two_j3*two_j4*two_j5*two_j6 == 0) then
         if (two_j1 == 0) then
            a = dble((two_j5 + 1)*(two_j3 + 1))
            ifa = two_j4 + two_j5 + two_j3
         else if (two_j2 == 0) then
            a = dble((two_j4 + 1)*(two_j3 + 1))
            ifa = two_j1 + two_j5 + two_j6
         else if (two_j3 == 0) then
            a = dble((two_j1 + 1)*(two_j4 + 1))
            ifa = two_j1 + two_j5 + two_j6
         else if (two_j4 == 0) then
            a = dble((two_j2 + 1)*(two_j3 + 1))
            ifa = two_j1 + two_j2 + two_j3
         else if (two_j5 == 0) then
            a = dble((two_j1 + 1)*(two_j3 + 1))
            ifa = two_j1 + two_j2 + two_j3
         else
            a = dble((two_j1 + 1)*(two_j2 + 1))
            ifa = two_j1 + two_j2 + two_j3
         endif
         six_j_value = one/dsqrt(a)
         if (mod(ifa,4) /= 0) six_j_value = -six_j_value
!
!     THE CASE 1/2
!
      else if (min0(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6) == 1) then
         if (two_j1 == 1) then
            call sixj5 (two_j5, two_j3, two_j4, two_j2, two_j6, 0, six_j_value)
         else if (two_j2 == 1) then
            call sixj5 (two_j1, two_j6, two_j5, two_j4, two_j3, 0, six_j_value)
         else if (two_j3 == 1) then
            call sixj5 (two_j1, two_j5, two_j6, two_j4, two_j2, 0, six_j_value)
         else if (two_j4 == 1) then
            call sixj5 (two_j2, two_j3, two_j1, two_j5, two_j6, 0, six_j_value)
         else if (two_j5 == 1) then
            call sixj5 (two_j1, two_j3, two_j2, two_j4, two_j6, 0, six_j_value)
         else
            call sixj5 (two_j1, two_j2, two_j3, two_j4, two_j5, 0, six_j_value)
         endif
!
!     THE CASE 1
!
      else if (min0(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6) == 2) then
         if (two_j1 == 2) then
            call sixj1 (two_j5, two_j3, two_j4, two_j2, two_j6, 0, six_j_value)
         else if (two_j2 == 2) then
            call sixj1 (two_j1, two_j6, two_j5, two_j4, two_j3, 0, six_j_value)
         else if (two_j3 == 2) then
            call sixj1 (two_j1, two_j5, two_j6, two_j4, two_j2, 0, six_j_value)
         else if (two_j4 == 2) then
            call sixj1 (two_j2, two_j3, two_j1, two_j5, two_j6, 0, six_j_value)
         else if (two_j5 == 2) then
            call sixj1 (two_j1, two_j3, two_j2, two_j4, two_j6, 0, six_j_value)
         else
            call sixj1 (two_j1, two_j2, two_j3, two_j4, two_j5, 0, six_j_value)
         endif
!
!     THE CASE 3/2
!
      else if (min0(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6) == 3) then
         if (two_j1 == 3) then
            call sixj35 (two_j5, two_j3, two_j4, two_j2, two_j6, 0, six_j_value)
         else if (two_j2 == 3) then
            call sixj35 (two_j1, two_j6, two_j5, two_j4, two_j3, 0, six_j_value)
         else if (two_j3 == 3) then
            call sixj35 (two_j1, two_j5, two_j6, two_j4, two_j2, 0, six_j_value)
         else if (two_j4 == 3) then
            call sixj35 (two_j2, two_j3, two_j1, two_j5, two_j6, 0, six_j_value)
         else if (two_j5 == 3) then
            call sixj35 (two_j1, two_j3, two_j2, two_j4, two_j6, 0, six_j_value)
         else
            call sixj35 (two_j1, two_j2, two_j3, two_j4, two_j5, 0, six_j_value)
         endif
!
!     THE CASE 2
!
      else if (min0(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6) == 4) then
         if (two_j1 == 4) then
            call sixj2 (two_j5, two_j3, two_j4, two_j2, two_j6, 0, six_j_value)
         else if (two_j2 == 4) then
            call sixj2 (two_j1, two_j6, two_j5, two_j4, two_j3, 0, six_j_value)
         else if (two_j3 == 4) then
            call sixj2 (two_j1, two_j5, two_j6, two_j4, two_j2, 0, six_j_value)
         else if (two_j4 == 4) then
            call sixj2 (two_j2, two_j3, two_j1, two_j5, two_j6, 0, six_j_value)
         else if (two_j5 == 4) then
            call sixj2 (two_j1, two_j3, two_j2, two_j4, two_j6, 0, six_j_value)
         else
            call sixj2 (two_j1, two_j2, two_j3, two_j4, two_j5, 0, six_j_value)
         endif
!
!     THE CASE 5/2
!
      else if (min0(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6) == 5) then
         call dracah (two_j1, two_j2, two_j5, two_j4, two_j3, two_j6, six_j_value)
         if (mod(two_j1 + two_j2 + two_j5 + two_j4,4) /= 0) six_j_value = -six_j_value
!
!     CASES 3
!
      else if (min0(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6) == 6) then
         if (two_j1 == 6) then
            call sixj3 (two_j5, two_j3, two_j4, two_j2, two_j6, 0, six_j_value)
         else if (two_j2 == 6) then
            call sixj3 (two_j1, two_j6, two_j5, two_j4, two_j3, 0, six_j_value)
         else if (two_j3 == 6) then
            call sixj3 (two_j1, two_j5, two_j6, two_j4, two_j2, 0, six_j_value)
         else if (two_j4 == 6) then
            call sixj3 (two_j2, two_j3, two_j1, two_j5, two_j6, 0, six_j_value)
         else if (two_j5 == 6) then
            call sixj3 (two_j1, two_j3, two_j2, two_j4, two_j6, 0, six_j_value)
         else
            call sixj3 (two_j1, two_j2, two_j3, two_j4, two_j5, 0, six_j_value)
         endif
!
!     THE CASE 7/2
!
      else if (min0(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6) == 7) then
         call dracah (two_j1, two_j2, two_j5, two_j4, two_j3, two_j6, six_j_value)
         if (mod(two_j1 + two_j2 + two_j5 + two_j4,4) /= 0) six_j_value = -six_j_value
!
!     CASES 4
!
      else if (min0(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6) == 8) then
         if (two_j1 == 8) then
            call sixj4 (two_j5, two_j3, two_j4, two_j2, two_j6, 0, six_j_value)
         else if (two_j2 == 8) then
            call sixj4 (two_j1, two_j6, two_j5, two_j4, two_j3, 0, six_j_value)
         else if (two_j3 == 8) then
            call sixj4 (two_j1, two_j5, two_j6, two_j4, two_j2, 0, six_j_value)
         else if (two_j4 == 8) then
            call sixj4 (two_j2, two_j3, two_j1, two_j5, two_j6, 0, six_j_value)
         else if (two_j5 == 8) then
            call sixj4 (two_j1, two_j3, two_j2, two_j4, two_j6, 0, six_j_value)
         else
            call sixj4 (two_j1, two_j2, two_j3, two_j4, two_j5, 0, six_j_value)
         endif
!
!     THE CASE 9/2
!
      else if (min0(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6) == 9) then
         call dracah (two_j1, two_j2, two_j5, two_j4, two_j3, two_j6, six_j_value)
         if (mod(two_j1 + two_j2 + two_j5 + two_j4,4) /= 0) six_j_value = -six_j_value
!
!     CALCULATED OTHER CASES
!
      else
       call dracah(two_j1,two_j2,two_j5,two_j4,two_j3,two_j6,six_j_value)
!         CALL GRACAH1 (I, J, M, L, K, N, SI)
!        CALL GRACAH(I,J,M,L,K,N,SI)
         if (mod(two_j1 + two_j2 + two_j5 + two_j4,4) /= 0) six_j_value = -six_j_value
      endif
      if (store_in_six_j_cache) small_six_j_cache(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6) = six_j_value
      return
      end subroutine sixj
end module sixj_i
