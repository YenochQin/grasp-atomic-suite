!*******************************************************************
!                                                                  *
      subroutine sixj1(i,j,k,l,m,itik,si)
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | I/2  J/2  K/2 |                                            *
!     | L/2  M/2   1  |               [B.M.X.  75].                *
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
      use cons_c,          only: zero, half, one, two, three
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ixjtik_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: i, j, k, l, m
      integer, intent(in)       :: itik
      real(real64), intent(out) :: si
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: ifa
      real(real64) :: as, aka, a, b, c
!-----------------------------------------------
      si = zero
      if (itik /= 0) then
!
!     CHESKED TRIANGULAR CONDITIONS
!
         if (ixjtik(i,j,k,l,m,2) == 0) return
      endif
      ifa = (i + j + k)/2
      as = dble(ifa)
      aka = one
      if (mod(ifa,2) /= 0) aka = -aka
      a = dble(k)
      b = dble(j)
      c = dble(i)
      if (i < m) then
         if (j < l) then
!              M > I,   L > J.
            si = aka*dsqrt((as + two)*(as + three)*(as - a + one)*(as - a + two&
               )/((b + one)*(b + two)*(b + three)*(c + one)*(c + two)*(c + &
               three)))
         else if (j == l) then
!              M > I,  L = J.
            si = (-aka)*dsqrt(two*(as + two)*(as - c)*(as - b + one)*(as - a + &
               one)/(b*(b + one)*(b + two)*(c + one)*(c + two)*(c + three)))
         else
!              M > I,  L < J.
            si = aka*dsqrt((as - c - one)*(as - c)*(as - b + one)*(as - b + two&
               )/((b - one)*b*(b + one)*(c + one)*(c + two)*(c + three)))
         endif
      else if (i == m) then
         if (j < l) then
!              M = L,  L > J.
            si = (-aka)*dsqrt((as + two)*(as - c + one)*(as - b)*(as - a + one)&
               *two/((b + one)*(b + two)*(b + three)*c*(c + one)*(c + two)))
         else if (j == l) then
!              M = I,  L = J.
            si = (-aka)*((b*b + c*c - a*a)*half + b + c - a)/dsqrt(b*(b + one)*&
               (b + two)*c*(c + one)*(c + two))
         else
!              M = I,  L < J.
            si = aka*dsqrt((as + one)*(as - c)*(as - b + one)*(as - a)*two/((b&
                - one)*b*(b + one)*c*(c + one)*(c + two)))
         endif
      else
         if (j < l) then
!              M < I,   L > J.
            si = aka*dsqrt((as - c + one)*(as - c + two)*(as - b - one)*(as - b&
               )/((b + one)*(b + two)*(b + three)*(c - one)*c*(c + one)))
         else if (j == l) then
!              M < I,   L = J.
            si = aka*dsqrt((as + one)*(as - c + one)*(as - b)*(as - a)*two/(b*(&
               b + one)*(b + two)*(c - one)*c*(c + one)))
         else
!              M < I,   L < J.
            si = aka*dsqrt(as*(as + one)*(as - a - one)*(as - a)/((b - one)*b*(&
               b + one)*(c - one)*c*(c + one)))
         endif
      endif
      return
      end subroutine sixj1
