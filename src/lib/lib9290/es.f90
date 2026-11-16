!***********************************************************************
!                                                                      *
      subroutine es(f, s2f, s3f)
!                                                                      *
!   Evaluate the sum of the series                                     *
!                                                                      *
!                       infinity      n              k                 *
!              S  (F) =   Sum     (-1)  exp (n*F) / n                  *
!               k        n = 0                                         *
!                                                                      *
!   for k = 2, 3 to machine precision. This is a utility subroutine,   *
!   called by SUBROUTINEs NUCPOT and NCHARG.                           *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 28 Sep 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:25   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(in) :: f
      real(real64), intent(out) :: s2f
      real(real64), intent(out) :: s3f
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: n
      real(real64) :: fase, en, obn, enf, term2, term3, s2last
!-----------------------------------------------
!
      n = 0
      s2f = 0.0_real64
      s3f = 0.0_real64
      fase = 1.0_real64
      n = n + 1
      en = dble(n)
      obn = 1.0_real64/en
      fase = -fase
      enf = exp(en*f)
      term2 = fase*enf*obn*obn
      term3 = term2*obn
      s2last = s2f
      s2f = s2f + term2
      s3f = s3f + term3
      do while(abs(s2f) /= abs(s2last))
         n = n + 1
         en = dble(n)
         obn = 1.0_real64/en
         fase = -fase
         enf = exp(en*f)
         term2 = fase*enf*obn*obn
         term3 = term2*obn
         s2last = s2f
         s2f = s2f + term2
         s3f = s3f + term3
      end do
!
      return
      end subroutine es
