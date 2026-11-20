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
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state,            only: pi,femtometer_to_bohr,speed_of_light_cm_per_s,hartree_to_wavenumber_cm
      use grasp_parameters,    only: max_radial_points, radial_workspace_capacity
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64) :: xvec(radial_workspace_capacity)
      real(real64) :: yvec(max_radial_points)
!GG      REAL(DOUBLE), DIMENSION(NNNP) :: XVEC, YVEC
      real(real64) :: p(5)
      real(real64) :: f(6)
      real(real64) :: par(2)
      real(real64) :: z, drms, rho, res
      integer      :: nrnuc, nparfit
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: x(radial_workspace_capacity)
      real(real64) :: y(max_radial_points), w(max_radial_points), dy(max_radial_points)
!GG      REAL(DOUBLE), DIMENSION(NNNP) :: X, Y, W, DY
      real(real64)  :: m(3, 3), mi(3, 3), c(3, 3), ci(3, 3)
      real(real64)    :: dr(4), dr2(4)
      real(real64)    :: rm(3), pm(3), b(3), df2(3)
      real(real64)    :: parf(2), parf2(2)
      real(real64) :: mdet, cdet, au2fm, norm, const, const2, a0
      real(real64) :: a1, a2, fdsum, fo90, r2,dr4,dr6, de, den, dden
      real(real64) :: dden2, dx2, dy2, x2, y2
      integer :: i, inuc, nmin, nmax, nr
!-----------------------------------------------
!
      au2fm = 1.0_real64/femtometer_to_bohr                ! Bohr radius in fm
      const = 1.0e-9_real64*speed_of_light_cm_per_s*hartree_to_wavenumber_cm*au2fm      ! To get electronic factors in GHz

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

      m(1,1) = 1.0_real64
      m(1,2) = x(nmin)**2.0_real64
      m(1,3) = x(nmin)**4.0_real64
      m(2,1) = 1.0_real64
      m(2,2) = x(nmin+2)**2.0_real64
      m(2,3) = x(nmin+2)**4.0_real64
      m(3,1) = 1.0_real64
      m(3,2) = x(nmin+4)**2.0_real64
      m(3,3) = x(nmin+4)**4.0_real64
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

      norm = -(y(nmax)-rho)/au2fm**3.0_real64

      do i=nmin,nmax
         x(i) = au2fm*x(i)
         y(i) = (y(i)-rho)/au2fm**3.0_real64
         y(i) = y(i)/norm
         w(i) = x(i)                        ! WEIGHTS SET TO R(I) TO COMPENSATE FOR EXPONENTIAL GRID DENSITY
      end do

!     THE DATA POINTS, SUBTRACTED SO THAT Y(I) = Y(I)-RHO, ARE FITTED TO POLYNOMIAL:
!     P(2)*X(I)**2 + P(3)*X(I)**4 + P(4)*X(I)**6
!     DETERMINE B_l AND C_{kl} MATRIX ELEMENTS
      b(:) = 0.0_real64
      c(:,:) = 0.0_real64
      do i=nmin,nmax
         b(1) = b(1) + y(i)*x(i)**2.0_real64*w(i)
         b(2) = b(2) + y(i)*x(i)**4.0_real64*w(i)
         b(3) = b(3) + y(i)*x(i)**6.0_real64*w(i)
         c(1,1) = c(1,1) + x(i)**4.0_real64*w(i)
         c(2,2) = c(2,2) + x(i)**8.0_real64*w(i)
         c(3,3) = c(3,3) + x(i)**12.0_real64*w(i)
         c(1,2) = c(1,2) + x(i)**6.0_real64*w(i)
         c(1,3) = c(1,3) + x(i)**8.0_real64*w(i)
         c(2,3) = c(2,3) + x(i)**10.0_real64*w(i)
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

      p(1) = rho/au2fm**3.0_real64
      p(2) = p(2)*norm
      p(3) = p(3)*norm
      p(4) = p(4)*norm

!     RESULTING ELECTRONIC FACTORS F(N) IN GHZ fm^{-2N}
      f(1) = 2.0_real64*pi/3.0_real64*z*p(1)*const
      f(2) = 2.0_real64*pi/10.0_real64*z*p(2)*const
      f(3) = 2.0_real64*pi/21.0_real64*z*p(3)*const
      f(4) = 2.0_real64*pi/36.0_real64*z*p(4)*const

!    DETERMINE AVERAGE POINT DISCREPANCY PARAMETER
      res = 0.0_real64
      do i=nmin,nmax
         y(i) = y(i)*au2fm**3.0_real64*norm+rho
         res = res + (y(i)-rho &
              -au2fm**3.0_real64*(p(2)*x(i)**2.0_real64+p(3)*x(i)**4.0_real64 &
              +p(4)*x(i)**6.0_real64))**2.0_real64
      end do

      res = sqrt(res/(nmax-nmin+1))
      res = res/rho*1000.0_real64    ! In per mille of RHO

!     BELOW WE APPROXIMATE THE IS ENERGY. DIFFERENCE IN RADIAL MOMENTS ARE CALUCLATED BETWEEN THE
!     CURRENT DISTRIBTUION (PARF) AND WITH A DISTRIBUTION WITH A R^2 1 FM LARGER. THE DIFFERNCE
!     IN RADIAL MOMENTS ARE THEN MULTIPLIED WITH THE ELECTRONIC FACTORS.

      parf(1) = par(1)*au2fm    ! 50% FALL OFF RADIUS C
      parf(2) = par(2)*au2fm    ! SKIN THICKNESS a. 10% TO 90% FALL OFF DISTANCE t = 4.0*ln(3)*a

      !DR2(1) = 1.0d0
      !PARF2(1) =  SQRT(1.d0+PARF(1)**2.d0)
      dr2(2) = 10.0_real64/7.0_real64*parf(1)**2.0_real64 + 30.0_real64/7.0_real64*(pi*parf(2))**2.0_real64
      dr2(3) = 5.0_real64/3.0_real64*parf(1)**4.0_real64 + &
           110.0_real64/9.0_real64*(pi*parf(2)*parf(1))**2.0_real64 + &
           239.0_real64/9.0_real64*(pi*parf(2))**4.0_real64
      dr2(4) = 20.0_real64/11.0_real64*parf(1)**6.0_real64 + &
           780.0_real64/33.0_real64*(pi*parf(2))**2.0_real64*parf(1)**4.0_real64 + &
           4100.0_real64/33.0_real64*(pi*parf(2))**4.0_real64*parf(1)**2.0_real64 + &
           8180.0_real64/33.0_real64*(pi*parf(2))**6.0_real64

      fdsum = f(1)
      do i=2,4
         fdsum = fdsum + f(i)*dr2(i)
      end do
      f(5) = fdsum              ! IS ENERGY IN GHZ / FM^2

      dr2(2) = 25.0_real64/21.0_real64
      dr2(3) = 25.0_real64/9.0_real64*parf(1)**2.0_real64 + &
           275.0_real64/27.0_real64*(pi*parf(2))**2.0_real64
      dr2(4) = 50.0_real64/11.0_real64*parf(1)**4.0_real64 +  &
           1300.0_real64/33.0_real64*(pi*parf(2)*parf(1))**2.0_real64 + &
           10250.0_real64/99.0_real64*(pi*parf(2))**4.0_real64

      fdsum = 0.0_real64
      do i=2,4
         fdsum = fdsum + f(i)*dr2(i)
      end do
      f(6) = fdsum              ! IS ENERGY IN GHZ / FM^2

      return

      end subroutine edensityfit
