!
!     ------------------------------------------------------------------
!     B N D I N V
!     ------------------------------------------------------------------
!
      subroutine bndinv(a, el, n, determ, epsil, itest, nsize)
!
!       DOUBLE PRECISION MATRIX INVERSION SUBROUTINE
!       FROM "DLYTAP".
!
!*      DOUBLE PRECISION E,F
!*      DOUBLE PRECISION A,EL,D,DSQRT,C,S,DETERP
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      integer , intent(out) :: itest
      integer , intent(in) :: nsize
      real(real64) , intent(out) :: determ
      real(real64) , intent(in) :: epsil
      real(real64) , intent(inout) :: a(nsize,nsize)
      real(real64) , intent(inout) :: el(nsize,nsize)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: isl2, k000fx, indsnl, i, j, n1, m, k, j1, i1, ks
      real(real64) :: d, c, s, deterp, f, e, epsilp, rat
!-----------------------------------------------
      if (n < 2) goto 140
      isl2 = 0
      k000fx = 2
      if (isl2 == 0) indsnl = 2
      if (isl2 == 1) indsnl = 1
!       CALL SLITET(2,INDSNL)
!       CALL OVERFL(K000FX)
!       CALL DVCHK(K000FX)
!
!       SET EL = IDENTITY MATRIX
      do i = 1, n
         el(i,:n) = 0.0_real64
         el(i,i) = 1.0_real64
      end do
!
!       TRIANGULARIZE A, FORM EL
!
      n1 = n - 1
      m = 2
      do j = 1, n1
         do i = m, n
            if (a(i,j) == 0.0_real64) cycle
            d = dsqrt(a(j,j)*a(j,j)+a(i,j)*a(i,j))
            c = a(j,j)/d
            s = a(i,j)/d
            do k = j, n
               d = c*a(j,k) + s*a(i,k)
               a(i,k) = c*a(i,k) - s*a(j,k)
               a(j,k) = d
            end do
            do k = 1, n
               d = c*el(j,k) + s*el(i,k)
               el(i,k) = c*el(i,k) - s*el(j,k)
               el(j,k) = d
            end do
         end do
         m = m + 1
      end do
!       CALL OVERFL(K000FX)
!       GO TO (140,51),K000FX
!
!       CALCULATE THE DETERMINANT
      deterp = a(1,1)
      do i = 2, n
         deterp = deterp*a(i,i)
      end do
      determ = deterp
!       CALL OVERFL(K000FX)
!       GO TO (140,520,520),K000FX
!
!       IS MATRIX SINGULAR
      f = a(1,1)
      e = a(1,1)
      do i = 2, n
         if (dabs(f) < dabs(a(i,i))) f = a(i,i)
         if (dabs(e) <= dabs(a(i,i))) cycle
         e = a(i,i)
      end do
      epsilp = epsil
      if (epsilp <= 0) epsilp = 1.0e-8_real64
      rat = e/f
      if (abs(rat) < epsilp) goto 130
!
!       INVERT TRIANGULAR MATRIX
      j = n
      do j1 = 1, n
!       CALL SLITE(2)
         i = j
         isl2 = 1
         do i1 = 1, j
!       CALL SLITET(2,K000FX)
            if (isl2 == 0) k000fx = 2
            if (isl2 == 1) then
               k000fx = 1
               isl2 = 0
            end if
            select case (k000fx)
            case default
               a(i,j) = 1.0_real64/a(i,i)
            case (2)
               ks = i + 1
               d = 0.0_real64
               d = sum(a(i,ks:j)*a(ks:j,j))
               a(i,j) = -d/a(i,i)
            end select
 1003       continue
            i = i - 1
         end do
         j = j - 1
      end do
!       CALL OVERFL(K000FX)
!       GO TO (140,103,103),K000FX

!103    CALL DVCHK(K000FX)
!       GO TO (140,105),K000FX
!
!       PREMULTIPLY EL BY INVERTED TRIANGULAR MATRIX
      m = 1
      do i = 1, n
         do j = 1, n
            d = 0.0_real64
            d = sum(a(i,m:n)*el(m:n,j))
            el(i,j) = d
         end do
         m = m + 1
      end do
!       CALL OVERFL(K000FX)
!       GO TO (140,123,123),K000FX
!
!       RECOPY EL TO A
      a(:n,:n) = el(:n,:n)
      itest = 0
!126    IF(INDSNL.EQ.1)CALL SLITE(2)
  126 continue
      if (indsnl == 1) isl2 = 1
      return
!
  130 continue
      itest = 1
      goto 126
  140 continue
      itest = -1
      goto 126
      end subroutine bndinv
