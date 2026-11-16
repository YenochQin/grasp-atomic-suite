!*******************************************************************
!                                                                  *
      subroutine sixj3(j, k, l, m, n, itik, si)
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
      use vast_kind_param, only: double
      use cons_c,          only: zero, one, two, three, four, seven
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ixjtik_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j, k, l, m, n
      integer, intent(in) :: itik
      real(double), intent(out) :: si
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: i1
      real(double) :: as, a, b, c, aka
!-----------------------------------------------
      si = zero
      if (itik /= 0) then
!
!     CHESKED TRIANGULAR CONDITIONS
!
         if (ixjtik(j,k,l,m,n,6) == 0) return
      endif
      i1 = (j + k + l)/2
      as = dble(i1)
      a = dble(l)
      b = dble(j)
      c = dble(k)
      aka = one
      if (mod(i1,2) /= 0) aka = -aka
      if (j - n == 6) then
! -3
         if (k - m == 6) then
!                        -3  -3
            si = aka*dsqrt((as + one)*as*(as - one)*(as - two)*(as - three)*(as&
                - four)/((b + one)*b*(b - one)*(b - two)*(b - three)*(b - four)&
               *(b - three - two)))
            si = si*dsqrt((as - a - three - two)*(as - a - four)*(as - a - &
               three)*(as - a - two)*(as - a - one)*(as - a)/((c + one)*c*(c - &
               one)*(c - two)*(c - three)*(c - four)*(c - three - two)))
         else if (m - k == 6) then
!                        3  -3
            si = aka*dsqrt((as - c)*(as - c - one)*(as - c - two)*(as - c - &
               three)*(as - c - four)*(as - c - three - two)/((b + one)*b*(b - &
               one)*(b - two)*(b - three)*(b - four)*(b - three - two)))
            si = si*dsqrt((as - b + three + three)*(as - b + three + two)*(as&
                - b + four)*(as - b + three)*(as - b + two)*(as - b + one)/((c&
                + four + three)*(c + three + three)*(c + three + two)*(c + four&
               )*(c + three)*(c + two)*(c + one)))
         else if (k - m == 4) then
!                       -2  -3
            si = aka*dsqrt(two*three*(as - c)*(as - b + one)*(as - three)*(as&
                - two)*(as - one)*as*(as + one)/((b - three - two)*(b - four)*(&
               b - three)*(b - two)*(b - one)*b*(b + one)))
            si = si*dsqrt((as - a - four)*(as - a - three)*(as - a - two)*(as&
                - a - one)*(as - a)/((c - four)*(c - three)*(c - two)*(c - one)&
               *c*(c + one)*(c + two)))
         else if (m - k == 4) then
!                        2  -3
            si = aka*dsqrt(two*three*(as + one)*(as - a)*(as - c - four)*(as - &
               c - three)*(as - c - two)*(as - c - one)*(as - c)/((c + three + &
               three)*(c + three + two)*(c + four)*(c + three)*(c + two)*(c + &
               one)*c))
            si = si*dsqrt((as - b + three + two)*(as - b + four)*(as - b + &
               three)*(as - b + two)*(as - b + one)/((b - three - two)*(b - &
               four)*(b - three)*(b - two)*(b - one)*b*(b + one)))
         else if (k - m == 2) then
!                       -1   -3
            si = aka*dsqrt(three*(three + two)*(as - c)*(as - c - one)*(as - b&
                + two)*(as - b + one)*(as - two)*(as - one)*as*(as + one)/((b&
                - three - two)*(b - four)*(b - three)*(b - two)*(b - one)*b*(b&
                + one)))
            si = si*dsqrt((as - a - three)*(as - a - two)*(as - a - one)*(as - &
               a)/((c - three)*(c - two)*(c - one)*c*(c + one)*(c + two)*(c + &
               three)))
         else if (m - k == 2) then
!                        1   -3
            si = aka*dsqrt(three*(three + two)*(as + one)*as*(as - a - one)*(as&
                - a)*(as - c - three)*(as - c - two)*(as - c - one)*(as - c)/((&
               b - three - two)*(b - four)*(b - three)*(b - two)*(b - one)*b*(b&
                + one)))
            si = si*dsqrt((as - b + four)*(as - b + three)*(as - b + two)*(as&
                - b + one)/((c + three + two)*(c + four)*(c + three)*(c + two)*&
               (c + one)*c*(c - one)))
         else if (m - k == 0) then
!                        0   -3
            si = aka*two*dsqrt((three + two)*(as + one)*as*(as - one)*(as - a&
                - two)*(as - a - one)*(as - a)/((b - three - two)*(b - four)*(b&
                - three)*(b - two)*(b - one)*b*(b + one)))
            si = si*dsqrt((as - c - two)*(as - c - one)*(as - c)*(as - b + &
               three)*(as - b + two)*(as - b + one)/((c + four)*(c + three)*(c&
                + two)*(c + one)*c*(c - one)*(c - two)))
         endif
      else if (n - j == 6) then
!  3
         if (k - m == 6) then
!                        -3  3
            si = aka*dsqrt((as - b)*(as - b - one)*(as - b - two)*(as - b - &
               three)*(as - b - four)*(as - b - three - two)/((c + one)*c*(c - &
               one)*(c - two)*(c - three)*(c - four)*(c - three - two)))
            si = si*dsqrt((as - c + three + three)*(as - c + three + two)*(as&
                - c + four)*(as - c + three)*(as - c + two)*(as - c + one)/((b&
                + four + three)*(b + three + three)*(b + three + two)*(b + four&
               )*(b + three)*(b + two)*(b + one)))
         else if (m - k == 6) then
!                        3   3
            si = aka*dsqrt((as - a + three + three)*(as - a + three + two)*(as&
                - a + four)*(as - a + three)*(as - a + two)*(as - a + one)/((b&
                + seven)*(b + three + three)*(b + three + two)*(b + four)*(b + &
               three)*(b + two)*(b + one)))
            si = si*dsqrt((as + seven)*(as + three + three)*(as + three + two)*&
               (as + four)*(as + three)*(as + two)/((c + seven)*(c + three + &
               three)*(c + three + two)*(c + four)*(c + three)*(c + two)*(c + &
               one)))
         else if (k - m == 4) then
!                       -2   3
            si = -aka*dsqrt(two*three*(as - a + one)*(as + two)*(as - b - four)&
               *(as - b - three)*(as - b - two)*(as - b - one)*(as - b)/((b + &
               seven)*(b + three + three)*(b + three + two)*(b + four)*(b + &
               three)*(b + two)*(b + one)))
            si = si*dsqrt((as - c + three + two)*(as - c + four)*(as - c + &
               three)*(as - c + two)*(as - c + one)/((c - four)*(c - three)*(c&
                - two)*(c - one)*c*(c + one)*(c + two)))
         else if (m - k == 4) then
!                        2   3
            si = -aka*dsqrt(two*three*(as - b)*(as - c + one)*(as - a + three&
                + two)*(as - a + four)*(as - a + three)*(as - a + two)*(as - a&
                + one)/((b + seven)*(b + three + three)*(b + three + two)*(b + &
               four)*(b + three)*(b + two)*(b + one)))
            si = si*dsqrt((as + three + three)*(as + three + two)*(as + four)*(&
               as + three)*(as + two)/((c + three + three)*(c + three + two)*(c&
                + four)*(c + three)*(c + two)*(c + one)*c))
         else if (k - m == 2) then
!                       -1   3
            si = aka*dsqrt(three*(three + two)*(as - a + one)*(as - a + two)*(&
               as + three)*(as + two)*(as - b - three)*(as - b - two)*(as - b&
                - one)*(as - b)/((b + seven)*(b + three + three)*(b + three + &
               two)*(b + four)*(b + three)*(b + two)*(b + one)))
            si = si*dsqrt((as - c + four)*(as - c + three)*(as - c + two)*(as&
                - c + one)/((c - three)*(c - two)*(c - one)*c*(c + one)*(c + &
               two)*(c + three)))
         else if (m - k == 2) then
!                        1   3
            si = aka*dsqrt(three*(three + two)*(as - b)*(as - b - one)*(as - c&
                + two)*(as - c + one)*(as - a + four)*(as - a + three)*(as - a&
                + two)*(as - a + one)/((b + seven)*(b + three + three)*(b + &
               three + two)*(b + four)*(b + three)*(b + two)*(b + one)))
            si = si*dsqrt((as + three + two)*(as + four)*(as + three)*(as + two&
               )/((c + three + two)*(c + four)*(c + three)*(c + two)*(c + one)*&
               c*(c - one)))
         else if (m - k == 0) then
!                        0   3
            si = -aka*two*dsqrt((three + two)*(as - b)*(as - b - one)*(as - b&
                - two)*(as - c + three)*(as - c + two)*(as - c + one)/((b + &
               seven)*(b + three + three)*(b + three + two)*(b + four)*(b + &
               three)*(b + two)*(b + one)))
            si = si*dsqrt((as - a + three)*(as - a + two)*(as - a + one)*(as + &
               four)*(as + three)*(as + two)/((c + four)*(c + three)*(c + two)*&
               (c + one)*c*(c - one)*(c - two)))
         endif
      else if (j - n == 4) then
! -2
         if (k - m == 6) then
!                       -3  -2
            si = aka*dsqrt(two*three*(as - b)*(as - c + one)*(as - three)*(as&
                - two)*(as - one)*as*(as + one)/((c - three - two)*(c - four)*(&
               c - three)*(c - two)*(c - one)*c*(c + one)))
            si = si*dsqrt((as - a - four)*(as - a - three)*(as - a - two)*(as&
                - a - one)*(as - a)/((b - four)*(b - three)*(b - two)*(b - one)&
               *b*(b + one)*(b + two)))
         else if (m - k == 6) then
!                       3   -2
            si = -aka*dsqrt(two*three*(as - a + one)*(as + two)*(as - c - four)&
               *(as - c - three)*(as - c - two)*(as - c - one)*(as - c)/((c + &
               seven)*(c + three + three)*(c + three + two)*(c + four)*(c + &
               three)*(c + two)*(c + one)))
            si = si*dsqrt((as - b + three + two)*(as - b + four)*(as - b + &
               three)*(as - b + two)*(as - b + one)/((b - four)*(b - three)*(b&
                - two)*(b - one)*b*(b + one)*(b + two)))
         else if (k - m == 4) then
!                      -2  -2
            si = aka*((two + three)*(as - c)*(as - b) - (as + two)*(as - a - &
               four))*dsqrt((as - two)*(as - one)*as*(as + one)/((b - four)*(b&
                - three)*(b - two)*(b - one)*b*(b + one)*(b + two)))
            si = si*dsqrt((as - a - three)*(as - a - two)*(as - a - one)*(as - &
               a)/((c - four)*(c - three)*(c - two)*(c - one)*c*(c + one)*(c + &
               two)))
         else if (m - k == 4) then
!                       2  -2
            si = -aka*((two + three)*(as + one)*(as - a + one) - (as - c + one)&
               *(as - b + three + two))*dsqrt((as - c - three)*(as - c - two)*(&
               as - c - one)*(as - c)/((b - four)*(b - three)*(b - two)*(b - &
               one)*b*(b + one)*(b + two)))
            si = si*dsqrt((as - b + four)*(as - b + three)*(as - b + two)*(as&
                - b + one)/((c + three + three)*(c + three + two)*(c + four)*(c&
                + three)*(c + two)*(c + one)*c))
         else if (k - m == 2) then
!                       -1  -2
            si = aka*(two*(as - c - one)*(as - b) - (as + two)*(as - a - three)&
               )*dsqrt(two*(three + two)*(as - c)*(as - b + one)*(as - one)*as*&
               (as + one)/((b - four)*(b - three)*(b - two)*(b - one)*b*(b + &
               one)*(b + two)))
            si = si*dsqrt((as - a - two)*(as - a - one)*(as - a)/((c - three)*(&
               c - two)*(c - one)*c*(c + one)*(c + two)*(c + three)))
         else if (m - k == 2) then
!                        1  -2
            si = -aka*(two*as*(as - a + one) - (as - c + one)*(as - b + four))*&
               dsqrt(two*(three + two)*(as + one)*(as - a)*(as - c - two)*(as&
                - c - one)*(as - c)/((b - four)*(b - three)*(b - two)*(b - one)&
               *b*(b + one)*(b + two)))
            si = si*dsqrt((as - b + three)*(as - b + two)*(as - b + one)/((c + &
               three + two)*(c + four)*(c + three)*(c + two)*(c + one)*c*(c - &
               one)))
         else if (m - k == 0) then
!                        0  -2
            si = -aka*((as - one)*(as - a + one) - (as - c + one)*(as - b + &
               three))*dsqrt(two*three*(three + two)*(as + one)*as*(as - a - &
               one)*(as - a)/((b - four)*(b - three)*(b - two)*(b - one)*b*(b&
                + one)*(b + two)))
            si = si*dsqrt((as - c - one)*(as - c)*(as - b + two)*(as - b + one)&
               /((c + four)*(c + three)*(c + two)*(c + one)*c*(c - one)*(c - &
               two)))
         endif
      else if (n - j == 4) then
!  2
         if (k - m == 6) then
!                        -3  2
            si = aka*dsqrt(two*three*(as + one)*(as - a)*(as - b - four)*(as - &
               b - three)*(as - b - two)*(as - b - one)*(as - b)/((b + three + &
               three)*(b + three + two)*(b + four)*(b + three)*(b + two)*(b + &
               one)*b))
            si = si*dsqrt((as - c + three + two)*(as - c + four)*(as - c + &
               three)*(as - c + two)*(as - c + one)/((c - three - two)*(c - &
               four)*(c - three)*(c - two)*(c - one)*c*(c + one)))
         else if (m - k == 6) then
!                        3  2
            si = -aka*dsqrt(two*three*(as - c)*(as - b + one)*(as - a + three&
                + two)*(as - a + four)*(as - a + three)*(as - a + two)*(as - a&
                + one)/((c + seven)*(c + three + three)*(c + three + two)*(c + &
               four)*(c + three)*(c + two)*(c + one)))
            si = si*dsqrt((as + three + three)*(as + three + two)*(as + four)*(&
               as + three)*(as + two)/((b + three + three)*(b + three + two)*(b&
                + four)*(b + three)*(b + two)*(b + one)*b))
         else if (k - m == 4) then
!                       -2  2
            si = -aka*((two + three)*(as + one)*(as - a + one) - (as - b + one)&
               *(as - c + three + two))*dsqrt((as - b - three)*(as - b - two)*(&
               as - b - one)*(as - b)/((c - four)*(c - three)*(c - two)*(c - &
               one)*c*(c + one)*(c + two)))
            si = si*dsqrt((as - c + four)*(as - c + three)*(as - c + two)*(as&
                - c + one)/((b + three + three)*(b + three + two)*(b + four)*(b&
                + three)*(b + two)*(b + one)*b))
         else if (m - k == 4) then
!                        2  2
            si = aka*((two + three)*(as - b)*(as - c) - (as - a)*(as + three + &
               three))*dsqrt((as - a + four)*(as - a + three)*(as - a + two)*(&
               as - a + one)/((b + three + three)*(b + three + two)*(b + four)*&
               (b + three)*(b + two)*(b + one)*b))
            si = si*dsqrt((as + three + two)*(as + four)*(as + three)*(as + two&
               )/((c + three + three)*(c + three + two)*(c + four)*(c + three)*&
               (c + two)*(c + one)*c))
         else if (k - m == 2) then
!                       -1  2
            si = aka*(two*(as - a + two)*(as + one) - (as - b + one)*(as - c + &
               four))*dsqrt(two*(three + two)*(as - a + one)*(as + two)*(as - b&
                - two)*(as - b - one)*(as - b)/((b + three + three)*(b + three&
                + two)*(b + four)*(b + three)*(b + two)*(b + one)*b))
            si = si*dsqrt((as - c + three)*(as - c + two)*(as - c + one)/((c - &
               three)*(c - two)*(c - one)*c*(c + one)*(c + two)*(c + three)))
         else if (m - k == 2) then
!                        1  2
            si = -aka*(two*(as - b - one)*(as - c) - (as - a)*(as + three + two&
               ))*dsqrt(two*(three + two)*(as - b)*(as - c + one)*(as - a + &
               three)*(as - a + two)*(as - a + one)/((b + three + three)*(b + &
               three + two)*(b + four)*(b + three)*(b + two)*(b + one)*b))
            si = si*dsqrt((as + four)*(as + three)*(as + two)/((c + three + two&
               )*(c + four)*(c + three)*(c + two)*(c + one)*c*(c - one)))
         else if (m - k == 0) then
!                        0  2
            si = aka*((as - b - two)*(as - c) - (as - a)*(as + four))*dsqrt(two&
               *three*(three + two)*(as - b)*(as - b - one)*(as - c + two)*(as&
                - c + one)/((b + three + three)*(b + three + two)*(b + four)*(b&
                + three)*(b + two)*(b + one)*b))
            si = si*dsqrt((as - a + two)*(as - a + one)*(as + three)*(as + two)&
               /((c + four)*(c + three)*(c + two)*(c + one)*c*(c - one)*(c - &
               two)))
         endif
      else if (j - n == 2) then
! - 1
         if (k - m == 6) then
!                       -3   -1
            si = aka*dsqrt(three*(three + two)*(as - b)*(as - b - one)*(as - c&
                + two)*(as - c + one)*(as - two)*(as - one)*as*(as + one)/((c&
                - three - two)*(c - four)*(c - three)*(c - two)*(c - one)*c*(c&
                + one)))
            si = si*dsqrt((as - a - three)*(as - a - two)*(as - a - one)*(as - &
               a)/((b - three)*(b - two)*(b - one)*b*(b + one)*(b + two)*(b + &
               three)))
         else if (m - k == 6) then
!                       3    -1
            si = aka*dsqrt(three*(three + two)*(as - a + one)*(as - a + two)*(&
               as + three)*(as + two)*(as - c - three)*(as - c - two)*(as - c&
                - one)*(as - c)/((c + seven)*(c + three + three)*(c + three + &
               two)*(c + four)*(c + three)*(c + two)*(c + one)))
            si = si*dsqrt((as - b + four)*(as - b + three)*(as - b + two)*(as&
                - b + one)/((b - three)*(b - two)*(b - one)*b*(b + one)*(b + &
               two)*(b + three)))
         else if (k - m == 4) then
!                       -2   -1
            si = aka*(two*(as - b - one)*(as - c) - (as + two)*(as - a - three)&
               )*dsqrt(two*(three + two)*(as - b)*(as - c + one)*(as - one)*as*&
               (as + one)/((c - four)*(c - three)*(c - two)*(c - one)*c*(c + &
               one)*(c + two)))
            si = si*dsqrt((as - a - two)*(as - a - one)*(as - a)/((b - three)*(&
               b - two)*(b - one)*b*(b + one)*(b + two)*(b + three)))
         else if (m - k == 4) then
!                       2    -1
            si = aka*(two*(as - a + two)*(as + one) - (as - c + one)*(as - b + &
               four))*dsqrt(two*(three + two)*(as - a + one)*(as + two)*(as - c&
                - two)*(as - c - one)*(as - c)/((c + three + three)*(c + three&
                + two)*(c + four)*(c + three)*(c + two)*(c + one)*c))
            si = si*dsqrt((as - b + three)*(as - b + two)*(as - b + one)/((b - &
               three)*(b - two)*(b - one)*b*(b + one)*(b + two)*(b + three)))
         else if (k - m == 2) then
!                        -1   -1
            si = aka*(two*three*(as - c)*(as - c - one)*(as - b - one)*(as - b)&
                - two*four*(as - c)*(as - b)*(as + two)*(as - a - two) + (as + &
               three)*(as + two)*(as - a - three)*(as - a - two))
            si = si*dsqrt(as*(as + one)*(as - a - one)*(as - a)/((b - three)*(b&
                - two)*(b - one)*b*(b + one)*(b + two)*(b + three)*(c - three)*&
               (c - two)*(c - one)*c*(c + one)*(c + two)*(c + three)))
         else if (m - k == 2) then
!                         1   -1
            si = aka*(two*three*(as + one)*as*(as - a + one)*(as - a + two) - &
               two*four*(as + one)*(as - a + one)*(as - c + one)*(as - b + &
               three) + (as - c + one)*(as - c + two)*(as - b + four)*(as - b&
                + three))
            si = si*dsqrt((as - c - one)*(as - c)*(as - b + two)*(as - b + one)&
               /((b - three)*(b - two)*(b - one)*b*(b + one)*(b + two)*(b + &
               three)*(c + three + two)*(c + four)*(c + three)*(c + two)*(c + &
               one)*c*(c - one)))
         else if (m - k == 0) then
!                         0   -1
            si = aka*two*(as*(as - one)*(as - a + one)*(as - a + two) - three*&
               as*(as - a + one)*(as - c + one)*(as - b + two) + (as - c + one)&
               *(as - c + two)*(as - b + three)*(as - b + two))
            si = si*dsqrt(three*(as + one)*(as - a)*(as - c)*(as - b + one)/((b&
                - three)*(b - two)*(b - one)*b*(b + one)*(b + two)*(b + three)*&
               (c + four)*(c + three)*(c + two)*(c + one)*c*(c - one)*(c - two)&
               ))
         endif
      else if (n - j == 2) then
! - 1
         if (k - m == 6) then
!                        -3   1
            si = aka*dsqrt(three*(three + two)*(as + one)*as*(as - a - one)*(as&
                - a)*(as - b - three)*(as - b - two)*(as - b - one)*(as - b)/((&
               c - three - two)*(c - four)*(c - three)*(c - two)*(c - one)*c*(c&
                + one)))
            si = si*dsqrt((as - c + four)*(as - c + three)*(as - c + two)*(as&
                - c + one)/((b + three + two)*(b + four)*(b + three)*(b + two)*&
               (b + one)*b*(b - one)))
         else if (m - k == 6) then
!                        3    1
            si = aka*dsqrt(three*(three + two)*(as - c)*(as - c - one)*(as - b&
                + two)*(as - b + one)*(as - a + four)*(as - a + three)*(as - a&
                + two)*(as - a + one)/((c + seven)*(c + three + three)*(c + &
               three + two)*(c + four)*(c + three)*(c + two)*(c + one)))
            si = si*dsqrt((as + three + two)*(as + four)*(as + three)*(as + two&
               )/((b + three + two)*(b + four)*(b + three)*(b + two)*(b + one)*&
               b*(b - one)))
         else if (k - m == 4) then
!                        -2   1
            si = -aka*(two*as*(as - a + one) - (as - b + one)*(as - c + four))*&
               dsqrt(two*(three + two)*(as + one)*(as - a)*(as - b - two)*(as&
                - b - one)*(as - b)/((c - four)*(c - three)*(c - two)*(c - one)&
               *c*(c + one)*(c + two)))
            si = si*dsqrt((as - c + three)*(as - c + two)*(as - c + one)/((b + &
               three + two)*(b + four)*(b + three)*(b + two)*(b + one)*b*(b - &
               one)))
         else if (m - k == 4) then
!                         2   1
            si = -aka*(two*(as - c - one)*(as - b) - (as - a)*(as + three + two&
               ))*dsqrt(two*(three + two)*(as - c)*(as - b + one)*(as - a + &
               three)*(as - a + two)*(as - a + one)/((c + three + three)*(c + &
               three + two)*(c + four)*(c + three)*(c + two)*(c + one)*c))
            si = si*dsqrt((as + four)*(as + three)*(as + two)/((b + three + two&
               )*(b + four)*(b + three)*(b + two)*(b + one)*b*(b - one)))
         else if (m - k == 2) then
!                         1   1
            si = aka*(two*three*(as - b)*(as - b - one)*(as - c)*(as - c - one)&
                - two*four*(as - b)*(as - c)*(as - a)*(as + four) + (as - a)*(&
               as - a - one)*(as + three + two)*(as + four))
            si = si*dsqrt((as - a + two)*(as - a + one)*(as + three)*(as + two)&
               /((b + three + two)*(b + four)*(b + three)*(b + two)*(b + one)*b&
               *(b - one)*(c + three + two)*(c + four)*(c + three)*(c + two)*(c&
                + one)*c*(c - one)))
         else if (k - m == 2) then
!                        -1   1
            si = aka*(two*three*(as + one)*as*(as - a + one)*(as - a + two) - &
               two*four*(as + one)*(as - a + one)*(as - b + one)*(as - c + &
               three) + (as - b + one)*(as - b + two)*(as - c + four)*(as - c&
                + three))
            si = si*dsqrt((as - b - one)*(as - b)*(as - c + two)*(as - c + one)&
               /((c - three)*(c - two)*(c - one)*c*(c + one)*(c + two)*(c + &
               three)*(b + three + two)*(b + four)*(b + three)*(b + two)*(b + &
               one)*b*(b - one)))
         else if (m - k == 0) then
!                         0   1
            si = -aka*two*((as - b - one)*(as - b - two)*(as - c)*(as - c - one&
               ) - three*(as - b - one)*(as - c)*(as - a)*(as + three) + (as - &
               a)*(as - a - one)*(as + four)*(as + three))
            si = si*dsqrt(three*(as - b)*(as - c + one)*(as - a + one)*(as + &
               two)/((b + three + two)*(b + four)*(b + three)*(b + two)*(b + &
               one)*b*(b - one)*(c + four)*(c + three)*(c + two)*(c + one)*c*(c&
                - one)*(c - two)))
         endif
      else if (j - n == 0) then
!  0
         if (k - m == 6) then
!                       -3    0
            si = aka*two*dsqrt((three + two)*(as + one)*as*(as - one)*(as - a&
                - two)*(as - a - one)*(as - a)/((c - three - two)*(c - four)*(c&
                - three)*(c - two)*(c - one)*c*(c + one)))
            si = si*dsqrt((as - b - two)*(as - b - one)*(as - b)*(as - c + &
               three)*(as - c + two)*(as - c + one)/((b + four)*(b + three)*(b&
                + two)*(b + one)*b*(b - one)*(b - two)))
         else if (m - k == 6) then
!                        3    0
            si = -aka*two*dsqrt((three + two)*(as - c)*(as - c - one)*(as - c&
                - two)*(as - b + three)*(as - b + two)*(as - b + one)/((c + &
               seven)*(c + three + three)*(c + three + two)*(c + four)*(c + &
               three)*(c + two)*(c + one)))
            si = si*dsqrt((as - a + three)*(as - a + two)*(as - a + one)*(as + &
               four)*(as + three)*(as + two)/((b + four)*(b + three)*(b + two)*&
               (b + one)*b*(b - one)*(b - two)))
         else if (k - m == 4) then
!                       -2    0
            si = -aka*((as - one)*(as - a + one) - (as - b + one)*(as - c + &
               three))*dsqrt(two*three*(three + two)*(as + one)*as*(as - a - &
               one)*(as - a)/((c - four)*(c - three)*(c - two)*(c - one)*c*(c&
                + one)*(c + two)))
            si = si*dsqrt((as - b - one)*(as - b)*(as - c + two)*(as - c + one)&
               /((b + four)*(b + three)*(b + two)*(b + one)*b*(b - one)*(b - &
               two)))
         else if (m - k == 4) then
!                        2    0
            si = aka*((as - c - two)*(as - b) - (as - a)*(as + four))*dsqrt(two&
               *three*(three + two)*(as - c)*(as - c - one)*(as - b + two)*(as&
                - b + one)/((c + three + three)*(c + three + two)*(c + four)*(c&
                + three)*(c + two)*(c + one)*c))
            si = si*dsqrt((as - a + two)*(as - a + one)*(as + three)*(as + two)&
               /((b + four)*(b + three)*(b + two)*(b + one)*b*(b - one)*(b - &
               two)))
         else if (k - m == 2) then
!                       -1    0
            si = aka*two*(as*(as - one)*(as - a + one)*(as - a + two) - three*&
               as*(as - a + one)*(as - b + one)*(as - c + two) + (as - b + one)&
               *(as - b + two)*(as - c + three)*(as - c + two))
            si = si*dsqrt(three*(as + one)*(as - a)*(as - b)*(as - c + one)/((c&
                - three)*(c - two)*(c - one)*c*(c + one)*(c + two)*(c + three)*&
               (b + four)*(b + three)*(b + two)*(b + one)*b*(b - one)*(b - two)&
               ))
         else if (m - k == 2) then
!                        1    0
            si = -aka*two*((as - c - one)*(as - c - two)*(as - b)*(as - b - one&
               ) - three*(as - c - one)*(as - b)*(as - a)*(as + three) + (as - &
               a)*(as - a - one)*(as + four)*(as + three))
            si = si*dsqrt(three*(as - c)*(as - b + one)*(as - a + one)*(as + &
               two)/((c + three + two)*(c + four)*(c + three)*(c + two)*(c + &
               one)*c*(c - one)*(b + four)*(b + three)*(b + two)*(b + one)*b*(b&
                - one)*(b - two)))
         else if (k - m == 0) then
!                        0    0
            si = aka*((as - b)*(as - b - one)*(as - b - two)*(as - c)*(as - c&
                - one)*(as - c - two) - three*three*(as - b)*(as - b - one)*(as&
                - c)*(as - c - one)*(as - a)*(as + two) + three*three*(as - b)*&
               (as - c)*(as - a)*(as - a - one)*(as + three)*(as + two) - (as&
                - a)*(as - a - one)*(as - a - two)*(as + four)*(as + three)*(as&
                + two))
            si = si/dsqrt((b + four)*(b + three)*(b + two)*(b + one)*b*(b - one&
               )*(b - two)*(c + four)*(c + three)*(c + two)*(c + one)*c*(c - &
               one)*(c - two))
         endif
      endif
      return
      end subroutine sixj3
