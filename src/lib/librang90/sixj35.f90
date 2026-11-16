!*******************************************************************
!                                                                  *
      subroutine sixj35(j,k,l,m,n,itik,si)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | J/2  K/2  L/2 |                                            *
!     | M/2  N/2  3/2 |             [B.M.X. 75]                    *
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
      use cons_c,          only: zero, one, two, three, four
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ixjtik_i
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
      integer :: i1
      real(real64) :: as, a, b, c, aka
!-----------------------------------------------
      si = zero
      if (itik /= 0) then
!
!     CHESKED TRIANGULAR CONDITIONS
!
         if (ixjtik(j,k,l,m,n,3) == 0) return
      endif
      i1 = (j + k + l)/2
      as = dble(i1)
      a = dble(l)
      b = dble(j)
      c = dble(k)
      aka = one
      if (mod(i1,2) /= 0) aka = -aka
      if (j - n == 3) then
! -3
         if (k - m == 3) then
!  I                      -3/2  -3/2
            si = aka*dsqrt((as - one)*as*(as + one)*(as - a - two)*(as - a - &
               one)*(as - a)/((b - two)*(b - one)*b*(b + one)*(c - two)*(c - &
               one)*c*(c + one)))
         else if (m - k == 3) then
!  IV  P(12)              3/2   -3/2
            si = aka*dsqrt((as - c - two)*(as - c - one)*(as - c)*(as - b + one&
               )*(as - b + two)*(as - b + three)/((c + 1)*(c + two)*(c + three)&
               *(c + four)*(b - two)*(b - one)*b*(b + one)))
         else if (k - m == 1) then
!  II  P(12)             -1/2   -3/2
            si = aka*dsqrt(three*as*(as + one)*(as - a - one)*(as - a)*(as - c)&
               *(as - b + one)/((c - one)*c*(c + one)*(c + two)*(b - two)*(b - &
               one)*b*(b + one)))
         else if (m - k == 1) then
!  III P(12)              1/2   -3/2
            si = aka*dsqrt(three*(as + one)*(as - a)*(as - c - one)*(as - c)*(&
               as - b + one)*(as - b + two)/(c*(c + one)*(c + two)*(c + three)*&
               (b - two)*(b - one)*b*(b + one)))
         endif
      else if (n - j == 3) then
!  3
         if (k - m == 3) then
!  IV                     -3/2   3/2
            si = aka*dsqrt((as - b - two)*(as - b - one)*(as - b)*(as - c + one&
               )*(as - c + two)*(as - c + three)/((b + one)*(b + two)*(b + &
               three)*(b + four)*(c - two)*(c - one)*c*(c + one)))
         else if (m - k == 3) then
!  2       pataisyta               3/2   3/2
            si = -aka*dsqrt((as + two)*(as + three)*(as + four)*(as - a + one)*&
               (as - a + two)*(as - a + three)/((b + one)*(b + two)*(b + three)&
               *(b + four)*(c + one)*(c + two)*(c + three)*(c + four)))
         else if (k - m == 1) then
!  1   P(12)   pataisytas          -1/2    3/2
            si = -aka*dsqrt(three*(as + two)*(as - a + one)*(as - c + one)*(as&
                - c + two)*(as - b - one)*(as - b)/((c - one)*c*(c + one)*(c + &
               two)*(b + one)*(b + two)*(b + three)*(b + four)))
         else if (m - k == 1) then
!  3  P(12)     taisyta           1/2    3/2
            si = aka*dsqrt(three*(as + two)*(as + three)*(as - a + one)*(as - a&
                + two)*(as - b)*(as - c + one)/(c*(c + one)*(c + two)*(c + &
               three)*(b + one)*(b + two)*(b + three)*(b + four)))
         endif
! -1
      else if (j - n == 1) then
         if (k - m == 3) then
!  II                   -3/2   -1/2
            si = aka*dsqrt((three*as*(as + one)*(as - a - one)*(as - a)*(as - b&
               )*(as - c + one))/((b - one)*b*(b + one)*(b + two)*(c - two)*(c&
                - one)*c*(c + one)))
         else if (m - k == 3) then
!  1                     3/2   -1/2
            si = -aka*dsqrt(three*(as + two)*(as - a + one)*(as - b + one)*(as&
                - b + two)*(as - c - one)*(as - c)/((b - one)*b*(b + one)*(b + &
               two)*(c + one)*(c + two)*(c + three)*(c + four)))
         else if (k - m == 1) then
!  V                    -1/2   -1/2
            si = aka*(two*(as - b)*(as - c) - (as + two)*(as - a - one))*dsqrt(&
               (as + one)*(as - a)/((b - one)*b*(b + one)*(b + two)*(c - one)*c&
               *(c + one)*(c + two)))
         else if (m - k == 1) then
!  VI P(12)              1/2   -1/2
            si = aka*((as - b + two)*(as - c + one) - two*(as - a + one)*(as + &
               one))*dsqrt((as - c)*(as - b + one)/(c*(c + one)*(c + two)*(c + &
               three)*(b - one)*b*(b + one)*(b + two)))
         endif
! 1
      else if (n - j == 1) then
         if (k - m == 3) then
!  III                  -3/2    1/2
            si = aka*dsqrt(three*(as + one)*(as - a)*(as - b - one)*(as - b)*(&
               as - c + one)*(as - c + two)/(b*(b + one)*(b + two)*(b + three)*&
               (c - two)*(c - one)*c*(c + one)))
         else if (m - k == 3) then
!  3              pataisyta       3/2    1/2
            si = aka*dsqrt(three*(as + two)*(as + three)*(as - a + one)*(as - a&
                + two)*(as - b + one)*(as - c)/(b*(b + one)*(b + two)*(b + &
               three)*(c + one)*(c + two)*(c + three)*(c + four)))
         else if (k - m == 1) then
!  VI                   -1/2    1/2
            si = aka*((as - c + two)*(as - b + one) - two*(as - a + one)*(as + &
               one))*dsqrt((as - b)*(as - c + one)/(b*(b + one)*(b + two)*(b + &
               three)*(c - one)*c*(c + one)*(c + two)))
         else if (m - k == 1) then
!  4      pataisyta               1/2    1/2
            si = -aka*(two*(as - b)*(as - c) - (as + three)*(as - a))*dsqrt((as&
                + two)*(as - a + one)/(b*(b + one)*(b + two)*(b + three)*c*(c&
                + one)*(c + two)*(c + three)))
         endif
      endif
      return
      end subroutine sixj35
