module sixj5_i
    implicit none (type, external)
    private
    public :: sixj5
contains
!*******************************************************************
!                                                                  *
      subroutine sixj5(j,k,l,m,n,itik,si)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | J/2  K/2  L/2 |                                            *
!     | M/2  N/2  1/2 |             [B.M.X. 75]                    *
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
      use cons_c,          only: zero, one, two
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ixjtik_i, only: ixjtik
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j, k, l, m, n
      integer, intent(in)       :: itik
      real(real64), intent(out) :: si
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: i1
      real(real64) :: as, a, b, c, aka
!-----------------------------------------------
      si = zero
      if (itik /= 0) then
!
!     CHESKED TRIANGULAR CONDITIONS
!
         if (ixjtik(j,k,l,m,n,1) == 0) return
      end if
      i1 = (j + k + l)/2
      as = dble(i1)
      a = dble(l)
      b = dble(k)
      c = dble(j)
      aka = one
      if (mod(i1,2) /= 0) aka = -aka
      if (k < m) then
         if (j < n) then
!              M > K,  J < N.
            si = -aka*dsqrt((as + two)*(as - a + one)/((b + one)*(b + two)*(c&
                + one)*(c + two)))
         else if (j > n) then
!              M > K,  J > N.
            si = aka*dsqrt((as - c + one)*(as - b)/((b + one)*(b + two)*c*(c + &
               one)))
         end if
      else if (k > m) then
         if (j < n) then
!             M < K,  J < N.
            si = aka*dsqrt((as - c)*(as - b + one)/(b*(b + one)*(c + one)*(c + &
               two)))
         else if (j > n) then
!             M < K,  J > N.
            si = aka*dsqrt((as + one)*(as - a)/(b*(b + one)*c*(c + one)))
         end if
      end if
      return
      end subroutine sixj5
end module sixj5_i
