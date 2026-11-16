!*******************************************************************
!                                                                  *
      subroutine sixj(i,j,k,l,m,n,itik,si)
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
      use vast_kind_param, only: double
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
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer             :: i, j, k, l, m, n
      integer, intent(in) :: itik
      real(double)        :: si
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ifa
      real(double), dimension(0:4,0:4,0:4,0:4,0:4,0:4) :: raca
      real(double) :: undef, a
      logical :: save
!     Persistent memoization is intentional (DATA already implies SAVE).
!     State remains per-thread when this source is compiled with OpenMP.
      save :: raca
!-----------------------------------------------
      data raca/ 15625*1.d-20/
      data undef/ 1.d-20/
!     The jj2lsj OpenMP object build gives each thread its own lazy cache.
!     Ordinary suite libraries compile this sentinel as a comment.
!$OMP THREADPRIVATE(RACA)
      si = zero
      if (itik /= 0) then
!
!     CHESKED TRIANGULAR CONDITIONS
!
         if (ixjtik(i,j,k,l,m,n) == 0) return
      endif
      save = .false.
      if (max0(i,j,k,l,m,n) <= 4) then
         si = raca(i,j,k,l,m,n)
         if (si == undef) then
            save = .true.
         else
            return
         endif
      endif
!
!     CALCULATED IN CASE WHEN ONE OF PERAMETERS EQUAL 0
!
      if (i*j*k*l*m*n == 0) then
         if (i == 0) then
            a = dble((m + 1)*(k + 1))
            ifa = l + m + k
         else if (j == 0) then
            a = dble((l + 1)*(k + 1))
            ifa = i + m + n
         else if (k == 0) then
            a = dble((i + 1)*(l + 1))
            ifa = i + m + n
         else if (l == 0) then
            a = dble((j + 1)*(k + 1))
            ifa = i + j + k
         else if (m == 0) then
            a = dble((i + 1)*(k + 1))
            ifa = i + j + k
         else
            a = dble((i + 1)*(j + 1))
            ifa = i + j + k
         endif
         si = one/dsqrt(a)
         if (mod(ifa,4) /= 0) si = -si
!
!     THE CASE 1/2
!
      else if (min0(i,j,k,l,m,n) == 1) then
         if (i == 1) then
            call sixj5 (m, k, l, j, n, 0, si)
         else if (j == 1) then
            call sixj5 (i, n, m, l, k, 0, si)
         else if (k == 1) then
            call sixj5 (i, m, n, l, j, 0, si)
         else if (l == 1) then
            call sixj5 (j, k, i, m, n, 0, si)
         else if (m == 1) then
            call sixj5 (i, k, j, l, n, 0, si)
         else
            call sixj5 (i, j, k, l, m, 0, si)
         endif
!
!     THE CASE 1
!
      else if (min0(i,j,k,l,m,n) == 2) then
         if (i == 2) then
            call sixj1 (m, k, l, j, n, 0, si)
         else if (j == 2) then
            call sixj1 (i, n, m, l, k, 0, si)
         else if (k == 2) then
            call sixj1 (i, m, n, l, j, 0, si)
         else if (l == 2) then
            call sixj1 (j, k, i, m, n, 0, si)
         else if (m == 2) then
            call sixj1 (i, k, j, l, n, 0, si)
         else
            call sixj1 (i, j, k, l, m, 0, si)
         endif
!
!     THE CASE 3/2
!
      else if (min0(i,j,k,l,m,n) == 3) then
         if (i == 3) then
            call sixj35 (m, k, l, j, n, 0, si)
         else if (j == 3) then
            call sixj35 (i, n, m, l, k, 0, si)
         else if (k == 3) then
            call sixj35 (i, m, n, l, j, 0, si)
         else if (l == 3) then
            call sixj35 (j, k, i, m, n, 0, si)
         else if (m == 3) then
            call sixj35 (i, k, j, l, n, 0, si)
         else
            call sixj35 (i, j, k, l, m, 0, si)
         endif
!
!     THE CASE 2
!
      else if (min0(i,j,k,l,m,n) == 4) then
         if (i == 4) then
            call sixj2 (m, k, l, j, n, 0, si)
         else if (j == 4) then
            call sixj2 (i, n, m, l, k, 0, si)
         else if (k == 4) then
            call sixj2 (i, m, n, l, j, 0, si)
         else if (l == 4) then
            call sixj2 (j, k, i, m, n, 0, si)
         else if (m == 4) then
            call sixj2 (i, k, j, l, n, 0, si)
         else
            call sixj2 (i, j, k, l, m, 0, si)
         endif
!
!     THE CASE 5/2
!
      else if (min0(i,j,k,l,m,n) == 5) then
         call dracah (i, j, m, l, k, n, si)
         if (mod(i + j + m + l,4) /= 0) si = -si
!
!     CASES 3
!
      else if (min0(i,j,k,l,m,n) == 6) then
         if (i == 6) then
            call sixj3 (m, k, l, j, n, 0, si)
         else if (j == 6) then
            call sixj3 (i, n, m, l, k, 0, si)
         else if (k == 6) then
            call sixj3 (i, m, n, l, j, 0, si)
         else if (l == 6) then
            call sixj3 (j, k, i, m, n, 0, si)
         else if (m == 6) then
            call sixj3 (i, k, j, l, n, 0, si)
         else
            call sixj3 (i, j, k, l, m, 0, si)
         endif
!
!     THE CASE 7/2
!
      else if (min0(i,j,k,l,m,n) == 7) then
         call dracah (i, j, m, l, k, n, si)
         if (mod(i + j + m + l,4) /= 0) si = -si
!
!     CASES 4
!
      else if (min0(i,j,k,l,m,n) == 8) then
         if (i == 8) then
            call sixj4 (m, k, l, j, n, 0, si)
         else if (j == 8) then
            call sixj4 (i, n, m, l, k, 0, si)
         else if (k == 8) then
            call sixj4 (i, m, n, l, j, 0, si)
         else if (l == 8) then
            call sixj4 (j, k, i, m, n, 0, si)
         else if (m == 8) then
            call sixj4 (i, k, j, l, n, 0, si)
         else
            call sixj4 (i, j, k, l, m, 0, si)
         endif
!
!     THE CASE 9/2
!
      else if (min0(i,j,k,l,m,n) == 9) then
         call dracah (i, j, m, l, k, n, si)
         if (mod(i + j + m + l,4) /= 0) si = -si
!
!     CALCULATED OTHER CASES
!
      else
       call dracah(i,j,m,l,k,n,si)
!         CALL GRACAH1 (I, J, M, L, K, N, SI)
!        CALL GRACAH(I,J,M,L,K,N,SI)
         if (mod(i + j + m + l,4) /= 0) si = -si
      endif
      if (save) raca(i,j,k,l,m,n) = si
      return
      end subroutine sixj
