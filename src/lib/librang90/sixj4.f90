!*******************************************************************
!                                                                  *
      subroutine sixj4(jc,je,jd,jb,jf,itik,si)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | JC/2  JE/2  JD/2 |                                         *
!     | JB/2  JF/2    4  |                                         *
!                                                                  *
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
      use cons_c,          only: zero, half, one, two, three, eps
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ixjtik_i
      use dracah_i
      use sixj2_i
      use sixj3_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jc, je, jd, jb, jf
      integer, intent(in) :: itik
      real(real64)        :: si
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: a, c, e, d, b, f, x1, x2, x3, s2, s3
!-----------------------------------------------
      si = zero
      if (itik /= 0) then
!
!     CHESKED TRIANGULAR CONDITIONS
!
         if (ixjtik(jc,je,jd,jb,jf,8) == 0) return
      endif
      if (ixjtik(jc,je,jd,jb,jf,6) == 0) then
         call dracah (jc, je, jf, jb, jd, 8, si)
         if (mod(jc + je + jf + jb,4) /= 0) si = -si
      else
         a = three
         c = dble(jc)*half
         e = dble(je)*half
         d = dble(jd)*half
         b = dble(jb)*half
         f = dble(jf)*half
         x1 = a*dsqrt((a+b+e+two)*(a-b+e+one)*(a+b-e+one)*((-         &
            a)+b+e)*(a+c+f+two)*(a-c+f+one)*(a+c-f+one)*(             &
            (-a)+c+f))
         x2 = (a+one)*dsqrt((a+b+e+one)*(a-b+e)*(a+b-e)*((-a)         &
             + b+e+one)*(a+c+f+one)*(a-c+f)*(a+c-f)*((-a)+c           &
             + f+one))
         x3 = (two*a+one)*(two*(a*(a+one)*d*(d+one)-b*(b+one)*c*(c+   &
            one)-e*(e+one)*f*(f+one))+(a*(a+one)-b*(b+one)-e*(e       &
             +one))*(a*(a+one)-c*(c+one)-f*(f+one)))
         if (dabs(x2) < eps) then
            s2 = zero
         else
            call sixj2 (jc, je, jd, jb, jf, 0, s2)
         endif
         if (dabs(x3) < eps) then
            s3 = zero
         else
            call sixj3 (jc, je, jd, jb, jf, 0, s3)
         endif
         si = (x3*s3 - x2*s2)/x1
      endif
      return
      end subroutine sixj4
