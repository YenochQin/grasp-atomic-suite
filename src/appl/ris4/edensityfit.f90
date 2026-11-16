!***********************************************************************
!                                                                      *
subroutine edensityfit(xvec,yvec,z,par,nrnuc,f,rho,res)
!SUBROUTINE EDENSITYFIT(XVEC,YVEC,Z,PAR,NPARFIT,                  &
!                                                DRMS,P,F,RHO,RES,NRNUC)

!     Fits polynomial b1 + b2r^2 b3r^3 + b4r^4 to (r,rho) electron     *
!     density data using least squares method                          *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use def_c,            only: pi,fmtoau,ccms,aucm
      use parameter_def,    only: nnnp, nnn1
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double), dimension(nnn1) :: xvec
      real(double), dimension(nnnp) :: yvec
!GG      REAL(DOUBLE), DIMENSION(NNNP) :: XVEC, YVEC
      real(double), dimension(5) :: p
      real(double), dimension(6) :: f
      real(double), dimension(2) :: par
      real(double) :: z, drms, rho, res
      integer      :: nrnuc, nparfit
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double), dimension(nnn1) :: x
      real(double), dimension(nnnp) :: y, w, dy
!GG      REAL(DOUBLE), DIMENSION(NNNP) :: X, Y, W, DY
      real(double), dimension(3,3)  :: m, mi, c, ci
      real(double), dimension(4)    :: dr, dr2
      real(double), dimension(3)    :: rm, pm, b, df2
      real(double), dimension(2)    :: parf, parf2
      real(double) :: mdet, cdet, au2fm, norm, const, const2, a0
      real(double) :: a1, a2, fdsum, fo90, r2,dr4,dr6, de, den, dden
      real(double) :: dden2, dx2, dy2, x2, y2
      integer :: i, inuc, nmin, nmax, nr
!-----------------------------------------------
!
      au2fm = 1.d0/fmtoau                ! Bohr radius in fm
      const = 1.d-9*ccms*aucm*au2fm      ! To get electronic factors in GHz

!     COPY ARRAYS
      x(:) = xvec(:)
      y(:) = yvec(:)

      ! DETERMINE FIRST POINT BEYOND FO90 IN GRID CALLED NR
      ! DETERMINE FIRST POINT WHERE GRID IS RELIABLE CALLED NMIN
      ! SINCE FIRST FEW DATA POINTS IN DENSITY ARE NOT RELIABLE DUE TO DIVISION WITH SMALL
      ! R^2 VALUES (STAGGERING IS SEEN), NMIN IS THE FIRST POINT TO BE USED IN THE SUBSEQUENT LEAST SQUARES FIT

      nmin = 10

!     SETS NMAX IN FIT TO NR
      nmax = nrnuc

!     DETERMINE RHO:
!     BELOW PM(1) + PM(2)*X(I)**2 + PM(3)*X(I)**4 IS FITTED TROUGH DATA POINTS
!     (X(NMIN),Y(NMIN)), (X(NMIN+2),Y(NMIN+2)), (X(NMIN+4),Y(NMIN+4)).
!     WE HAVE: RHO(0) = PM(1)

      rm(1) = y(nmin)
      rm(2) = y(nmin+2)
      rm(3) = y(nmin+4)

      m(1,1) = 1.0d0
      m(1,2) = x(nmin)**2.0d0
      m(1,3) = x(nmin)**4.0d0
      m(2,1) = 1.0d0
      m(2,2) = x(nmin+2)**2.0d0
      m(2,3) = x(nmin+2)**4.0d0
      m(3,1) = 1.d0
      m(3,2) = x(nmin+4)**2.0d0
      m(3,3) = x(nmin+4)**4.0d0
      mdet = m(1,1)*m(2,2)*m(3,3)-m(1,1)*m(2,3)*m(3,2)-                &
           m(1,2)*m(2,1)*m(3,3)+m(1,2)*m(2,3)*m(3,1)+                  &
           m(1,3)*m(2,1)*m(3,2)-m(1,3)*m(2,2)*m(3,1)

! Determine inverse of C matrix
      mi(1,1) = (m(2,2)*m(3,3)-m(2,3)*m(3,2))/mdet
      mi(1,2) = (m(1,3)*m(3,2)-m(1,2)*m(3,3))/mdet
      mi(1,3) = (m(1,2)*m(2,3)-m(1,3)*m(2,2))/mdet

      mi(2,1) = (m(2,3)*m(3,1)-m(2,1)*m(3,3))/mdet
      mi(2,2) = (m(1,1)*m(3,3)-m(1,3)*m(3,1))/mdet
      mi(2,3) = (m(1,3)*m(2,1)-m(1,1)*m(2,3))/mdet

      mi(3,1) = (m(2,1)*m(3,2)-m(2,2)*m(3,1))/mdet
      mi(3,2) = (m(1,2)*m(3,1)-m(1,1)*m(3,2))/mdet
      mi(3,3) = (m(1,1)*m(2,2)-m(1,2)*m(2,1))/mdet

      ! Determine parameters
      pm(1) = mi(1,1)*rm(1)+mi(1,2)*rm(2)+mi(1,3)*rm(3)
      pm(2) = mi(2,1)*rm(1)+mi(2,2)*rm(2)+mi(2,3)*rm(3)
      pm(3) = mi(3,1)*rm(1)+mi(3,2)*rm(2)+mi(3,3)*rm(3)

      ! Finallay RHO(0) in au^{-3} is determined
      rho = pm(1)

      ! START LEAST SQUARE FIT PROCEDURE FOR DATA POINTS
      ! (X(NMIN),Y(NMIN)), (X(NMIN+1),Y(NMIN+1)),...,(X(N),Y(N))

      norm = -(y(nmax)-rho)/au2fm**3.0d0

      do i=nmin,nmax
         x(i) = au2fm*x(i)
         y(i) = (y(i)-rho)/au2fm**3.0d0
         y(i) = y(i)/norm
         w(i) = x(i)                        ! WEIGHTS SET TO R(I) TO COMPENSATE FOR EXPONENTIAL GRID DENSITY
      end do

!     THE DATA POINTS, SUBTRACTED SO THAT Y(I) = Y(I)-RHO, ARE FITTED TO POLYNOMIAL:
!     P(2)*X(I)**2 + P(3)*X(I)**4 + P(4)*X(I)**6
!     DETERMINE B_l AND C_{kl} MATRIX ELEMENTS
      b(:) = 0.0d0
      c(:,:) = 0.0d0
      do i=nmin,nmax
         b(1) = b(1) + y(i)*x(i)**2.0d0*w(i)
         b(2) = b(2) + y(i)*x(i)**4.0d0*w(i)
         b(3) = b(3) + y(i)*x(i)**6.0d0*w(i)
         c(1,1) = c(1,1) + x(i)**4.0d0*w(i)
         c(2,2) = c(2,2) + x(i)**8.0d0*w(i)
         c(3,3) = c(3,3) + x(i)**12.0d0*w(i)
         c(1,2) = c(1,2) + x(i)**6.0d0*w(i)
         c(1,3) = c(1,3) + x(i)**8.0d0*w(i)
         c(2,3) = c(2,3) + x(i)**10.0d0*w(i)
      end do
      c(2,1) = c(1,2)
      c(3,1) = c(1,3)
      c(3,2) = c(2,3)

!     COMPUTE DETERMINANT
      cdet = c(1,1)*c(2,2)*c(3,3)-c(1,1)*c(2,3)*c(3,2)- &
           c(1,2)*c(2,1)*c(3,3)+c(1,2)*c(2,3)*c(3,1)+ &
           c(1,3)*c(2,1)*c(3,2)-c(1,3)*c(2,2)*c(3,1)

!     COMPUTE INVERSE OF C MATRIX
      ci(1,1) = (c(2,2)*c(3,3)-c(2,3)*c(3,2))/cdet
      ci(1,2) = (c(1,3)*c(3,2)-c(1,2)*c(3,3))/cdet
      ci(1,3) = (c(1,2)*c(2,3)-c(1,3)*c(2,2))/cdet

      ci(2,1) = (c(2,3)*c(3,1)-c(2,1)*c(3,3))/cdet
      ci(2,2) = (c(1,1)*c(3,3)-c(1,3)*c(3,1))/cdet
      ci(2,3) = (c(1,3)*c(2,1)-c(1,1)*c(2,3))/cdet

      ci(3,1) = (c(2,1)*c(3,2)-c(2,2)*c(3,1))/cdet
      ci(3,2) = (c(1,2)*c(3,1)-c(1,1)*c(3,2))/cdet
      ci(3,3) = (c(1,1)*c(2,2)-c(1,2)*c(2,1))/cdet

!     DETERMINE FITTING PARAMETERS
      p(2) = ci(1,1)*b(1)+ci(1,2)*b(2)+ci(1,3)*b(3)
      p(3) = ci(2,1)*b(1)+ci(2,2)*b(2)+ci(2,3)*b(3)
      p(4) = ci(3,1)*b(1)+ci(3,2)*b(2)+ci(3,3)*b(3)

      p(1) = rho/au2fm**3.0d0
      p(2) = p(2)*norm
      p(3) = p(3)*norm
      p(4) = p(4)*norm

!     RESULTING ELECTRONIC FACTORS F(N) IN GHZ fm^{-2N}
      f(1) = 2.0d0*pi/3.0d0*z*p(1)*const
      f(2) = 2.0d0*pi/10.0d0*z*p(2)*const
      f(3) = 2.0d0*pi/21.0d0*z*p(3)*const
      f(4) = 2.0d0*pi/36.0d0*z*p(4)*const

!    DETERMINE AVERAGE POINT DISCREPANCY PARAMETER
      res = 0.0d0
      do i=nmin,nmax
         y(i) = y(i)*au2fm**3.0d0*norm+rho
         res = res + (y(i)-rho &
              -au2fm**3.0d0*(p(2)*x(i)**2.0d0+p(3)*x(i)**4.0d0 &
              +p(4)*x(i)**6.0d0))**2.0d0
      end do

      res = sqrt(res/(nmax-nmin+1))
      res = res/rho*1000.0d0    ! In per mille of RHO

!     BELOW WE APPROXIMATE THE IS ENERGY. DIFFERENCE IN RADIAL MOMENTS ARE CALUCLATED BETWEEN THE
!     CURRENT DISTRIBTUION (PARF) AND WITH A DISTRIBUTION WITH A R^2 1 FM LARGER. THE DIFFERNCE
!     IN RADIAL MOMENTS ARE THEN MULTIPLIED WITH THE ELECTRONIC FACTORS.

      parf(1) = par(1)*au2fm    ! 50% FALL OFF RADIUS C
      parf(2) = par(2)*au2fm    ! SKIN THICKNESS a. 10% TO 90% FALL OFF DISTANCE t = 4.0*ln(3)*a

      !DR2(1) = 1.0d0
      !PARF2(1) =  SQRT(1.d0+PARF(1)**2.d0)
      dr2(2) = 10.d0/7.d0*parf(1)**2.d0 + 30.d0/7.d0*(pi*parf(2))**2.d0
      dr2(3) = 5.d0/3.d0*parf(1)**4.d0 + &
           110.d0/9.d0*(pi*parf(2)*parf(1))**2.d0 + &
           239.d0/9.d0*(pi*parf(2))**4.d0
      dr2(4) = 20.d0/11.d0*parf(1)**6.d0 + &
           780.d0/33.d0*(pi*parf(2))**2.d0*parf(1)**4.d0 + &
           4100.d0/33.d0*(pi*parf(2))**4.d0*parf(1)**2.d0 + &
           8180.d0/33.d0*(pi*parf(2))**6.d0

      fdsum = f(1)
      do i=2,4
         fdsum = fdsum + f(i)*dr2(i)
      end do
      f(5) = fdsum              ! IS ENERGY IN GHZ / FM^2

      dr2(2) = 25.d0/21.d0
      dr2(3) = 25.d0/9.d0*parf(1)**2.d0 + &
           275.d0/27.d0*(pi*parf(2))**2.d0
      dr2(4) = 50.d0/11.d0*parf(1)**4.d0 +  &
           1300.d0/33.d0*(pi*parf(2)*parf(1))**2.d0 + &
           10250.d0/99.d0*(pi*parf(2))**4.d0

      fdsum = 0.d0
      do i=2,4
         fdsum = fdsum + f(i)*dr2(i)
      end do
      f(6) = fdsum              ! IS ENERGY IN GHZ / FM^2

      return

      end subroutine edensityfit
