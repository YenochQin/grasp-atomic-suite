module sixj2_i
    implicit none (type, external)
    private
    public :: sixj2
contains
!*******************************************************************
!                                                                  *
      subroutine sixj2(j,k,l,m,n,itik,si)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | J/2  K/2  L/2 |                                            *
!     | M/2  N/2   2  |             [B.M.X. 75]                    *
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
      use cons_c,          only: zero, half, one, two, three, four
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
      real(real64) :: as, a, b, c, aka, x1, x2, x3
!-----------------------------------------------
      si = zero
      if (itik /= 0) then
!
!     CHESKED TRIANGULAR CONDITIONS
!
         if (ixjtik(j,k,l,m,n,4) == 0) return
      end if
      i1 = (j + k + l)/2
      as = dble(i1)
      a = dble(l)
      b = dble(j)
      c = dble(k)
      aka = one
      if (mod(i1,2) /= 0) aka = -aka
      if (j - n == 4) then
! -2
         if (k - m == 4) then
!  I                      -2  -2
            si = aka*dsqrt((as - two)*(as - one)*as*(as + one)/((b - three)*(b&
                - two)*(b - one)*b*(b + one)))
            si = si*dsqrt((as - a - three)*(as - a - two)*(as - a - one)*(as - &
               a)/((c - three)*(c - two)*(c - one)*c*(c + one)))
         else if (m - k == 4) then
!  V   P(12)               2  -2
            si = aka*dsqrt((as - c - three)*(as - c - two)*(as - c - one)*(as&
                - c)/((c + one)*(c + two)*(c + three)*(c + four)*(c + two + &
               three)))
            si = si*dsqrt((as - b + one)*(as - b + two)*(as - b + three)*(as - &
               b + four)/((b - three)*(b - two)*(b - one)*b*(b + one)))
         else if (k - m == 2) then
!  II   P(12)              -1   -2
            si = aka*two*dsqrt((as - one)*as*(as + one)/((c - two)*(c - one)*c*&
               (c + one)*(c + two)))
            si = si*dsqrt((as - a - two)*(as - a - one)*(as - a)*(as - c)*(as&
                - b + one)/((b - three)*(b - two)*(b - one)*b*(b + one)))
         else if (m - k == 2) then
!  IV   P(12)               1   -2
            si = aka*two*dsqrt((as + one)*(as - a)*(as - c - two)*(as - c - one&
               )*(as - c)/(c*(c + one)*(c + two)*(c + three)*(c + four)))
            si = si*dsqrt((as - b + one)*(as - b + two)*(as - b + three)/((b - &
               three)*(b - two)*(b - one)*b*(b + one)))
         else if (k - m == 0) then
!  III  P(12)               0   -2
            si = aka*dsqrt(two*three*as*(as + one)*(as - a - one)*(as - a)/((c&
                - one)*c*(c + one)*(c + two)*(c + three)))
            si = si*dsqrt((as - c - one)*(as - c)*(as - b + one)*(as - b + two)&
               /((b - three)*(b - two)*(b - one)*b*(b + one)))
         end if
!  2
      else if (n - j == 4) then
         if (k - m == 4) then
!  V                      -2   2
            si = aka*dsqrt((as - b - three)*(as - b - two)*(as - b - one)*(as&
                - b)/((b + one)*(b + two)*(b + three)*(b + four)*(b + two + &
               three)))
            si = si*dsqrt((as - c + one)*(as - c + two)*(as - c + three)*(as - &
               c + four)/((c - three)*(c - two)*(c - one)*c*(c + one)))
         else if (m - k == 4) then
!  1                       2   2
            si = aka*dsqrt((as - a + four)*(as - a + three)*(as - a + two)*(as&
                - a + one)/((b + three + two)*(b + four)*(b + three)*(b + two)*&
               (b + one)))
            si = si*dsqrt((as + three + two)*(as + four)*(as + three)*(as + two&
               )/((c + three + two)*(c + four)*(c + three)*(c + two)*(c + one))&
               )
         else if (k - m == 2) then
!  3                      -1   2
            si = -aka*dsqrt((as - a + one)*(as + two)*(as - b - two)*(as - b - &
               one)*(as - b)/((b + two + three)*(b + four)*(b + three)*(b + two&
               )*(b + one)))
            si = si*two*dsqrt((as - c + three)*(as - c + two)*(as - c + one)/((&
               c - two)*(c - one)*c*(c + one)*(c + two)))
         else if (m - k == 2) then
!  2                       1   2
            si = -aka*dsqrt((as - b)*(as - c + one)*(as - a + three)*(as - a + &
               two)*(as - a + one)/((b + three + two)*(b + four)*(b + three)*(b&
                + two)*(b + one)))
            si = si*two*dsqrt((as + four)*(as + three)*(as + two)/((c + four)*(&
               c + three)*(c + two)*(c + one)*c))
         else if (k - m == 0) then
!  5                        0   2
            si = aka*dsqrt(three*two*(as - b)*(as - b - one)*(as - c + two)*(as&
                - c + one)/((b + three + two)*(b + four)*(b + three)*(b + two)*&
               (b + one)))
            si = si*dsqrt((as - a + two)*(as - a + one)*(as + three)*(as + two)&
               /((c + three)*(c + two)*(c + one)*c*(c - one)))
         end if
      else if (j - n == 2) then
! -1
         if (k - m == 4) then
!  II   P(12)              -2  -1
            si = aka*two*dsqrt((as - one)*as*(as + one)/((b - two)*(b - one)*b*&
               (b + one)*(b + two)))
            si = si*dsqrt((as - a - two)*(as - a - one)*(as - a)*(as - b)*(as&
                - c + one)/((c - three)*(c - two)*(c - one)*c*(c + one)))
         else if (m - k == 4) then
!  3   P(12)                2  -1
            si = -aka*dsqrt((as - a + one)*(as + two)*(as - c - two)*(as - c - &
               one)*(as - c)/((c + two + three)*(c + four)*(c + three)*(c + two&
               )*(c + one)))
            si = si*two*dsqrt((as - b + three)*(as - b + two)*(as - b + one)/((&
               b - two)*(b - one)*b*(b + one)*(b + two)))
         else if (k - m == 2) then
!  VI                       -1  -1
            si = aka*((a + b)*(a - b + two) - (c - two)*(c - b + two))/dsqrt((b&
                - two)*(b - one)*b*(b + one)*(b + two))
            si = si*dsqrt(as*(as + one)*(as - a - one)*(as - a)/((c - two)*(c&
                - one)*c*(c + one)*(c + two)))
         else if (m - k == 2) then
!  VIII P(12)              1  -1
            si = aka*((a + c + four)*(a - c - two) - (b - two)*(b + c + four))/&
               dsqrt(c*(c + one)*(c + two)*(c + three)*(c + four))
            si = si*dsqrt((as - c - one)*(as - c)*(as - b + one)*(as - b + two)&
               /((b - two)*(b - one)*b*(b + one)*(b + two)))
         else if (k - m == 0) then
!  VII  P(12)              0  -1
            si = aka*half*((a + c + two)*(a - c) - b*b + four)/dsqrt((c - one)*&
               c*(c + one)*(c + two)*(c + three))
            si = si*dsqrt(three*two*(as + one)*(as - a)*(as - c)*(as - b + one)&
               /((b - two)*(b - one)*b*(b + one)*(b + two)))
         end if
      else if (n - j == 2) then
!  1
         if (k - m == 4) then
!  IV                     -2   1
            si = aka*two*dsqrt((as + one)*(as - a)*(as - b - two)*(as - b - one&
               )*(as - b)/(b*(b + one)*(b + two)*(b + three)*(b + four)))
            si = si*dsqrt((as - c + one)*(as - c + two)*(as - c + three)/((c - &
               three)*(c - two)*(c - one)*c*(c + one)))
         else if (m - k == 4) then
!  2 P(12)                 2   1
            si = -aka*dsqrt((as - c)*(as - b + one)*(as - a + three)*(as - a + &
               two)*(as - a + one)/((c + three + two)*(c + four)*(c + three)*(c&
                + two)*(c + one)))
            si = si*two*dsqrt((as + four)*(as + three)*(as + two)/((b + four)*(&
               b + three)*(b + two)*(b + one)*b))
         else if (k - m == 2) then
!  VIII                   -1   1
            si = aka*((a + b + four)*(a - b - two) - (c - two)*(b + c + four))/&
               dsqrt(b*(b + one)*(b + two)*(b + three)*(b + four))
            si = si*dsqrt((as - b - one)*(as - b)*(as - c + one)*(as - c + two)&
               /((c - two)*(c - one)*c*(c + one)*(c + two)))
         else if (m - k == 2) then
!  4                       1   1
            si = aka*(three*(as - b)*(as - c) - (as - a)*(as + four))/dsqrt((b&
                + four)*(b + three)*(b + two)*(b + one)*b)
            si = si*dsqrt((as - a + two)*(as - a + one)*(as + three)*(as + two)&
               /((c + four)*(c + three)*(c + two)*(c + one)*c))
         else if (k - m == 0) then
!  6 P(12)                 0   1
            si = -aka*((as - b - one)*(as - c) - (as - a)*(as + three))/dsqrt((&
               b + four)*(b + three)*(b + two)*(b + one)*b)
            si = si*dsqrt(three*two*(as - b)*(as - c + one)*(as - a + one)*(as&
                + two)/((c + three)*(c + two)*(c + one)*c*(c - one)))
         end if
      else if (n - j == 0) then
! 0
         if (k - m == 4) then
!  III                     -2   0
            si = aka*dsqrt(three*two*as*(as + one)*(as - a - one)*(as - a)/((b&
                - one)*b*(b + one)*(b + two)*(b + three)))
            si = si*dsqrt((as - b - one)*(as - b)*(as - c + one)*(as - c + two)&
               /((c - three)*(c - two)*(c - one)*c*(c + one)))
         else if (m - k == 4) then
!  5                        2   0
            si = aka*dsqrt(three*two*(as - c)*(as - c - one)*(as - b + two)*(as&
                - b + one)/((c + three + two)*(c + four)*(c + three)*(c + two)*&
               (c + one)))
            si = si*dsqrt((as - a + two)*(as - a + one)*(as + three)*(as + two)&
               /((b + three)*(b + two)*(b + one)*b*(b - one)))
         else if (k - m == 2) then
!  VII                     -1   0
            si = aka*half*((a + b + two)*(a - b) - c*c + four)/dsqrt((b - one)*&
               b*(b + one)*(b + two)*(b + three))
            si = si*dsqrt(three*two*(as + one)*(as - a)*(as - b)*(as - c + one)&
               /((c - two)*(c - one)*c*(c + one)*(c + two)))
         else if (m - k == 2) then
!  6                        1   0
            si = -aka*((as - c - one)*(as - b) - (as - a)*(as + three))/dsqrt((&
               c + four)*(c + three)*(c + two)*(c + one)*c)
            si = si*dsqrt(three*two*(as - c)*(as - b + one)*(as - a + one)*(as&
                + two)/((b + three)*(b + two)*(b + one)*b*(b - one)))
         else if (k - m == 0) then
!  IX                       0   0
            x1 = (as - b)*(as - b - one)*(as - c)*(as - c - one)
            x2 = four*(as - b)*(as - c)*(as - a)*(as + two)
            x3 = (as - a)*(as - a - one)*(as + three)*(as + two)
            si = aka*(x1 - x2 + x3)/dsqrt((b - one)*b*(b + one)*(b + two)*(b + &
               three)*(c - one)*c*(c + one)*(c + two)*(c + three))
         end if
      end if
      return
      end subroutine sixj2
end module sixj2_i
