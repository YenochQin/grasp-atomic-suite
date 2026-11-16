module dcbsrw_i
    implicit none (type, external)
    private
    public :: dcbsrw
contains
!***********************************************************************
!                                                                      *
      subroutine dcbsrw(n, kappa, z, e, rg0, rg, rf, mtp)
!                                                                      *
!   This subroutine computes the  Dirac-Coulomb  bound-state orbital   *
!   radial wavefunction.   Equations (13.5) and (13.5') of  Akhiezer   *
!   and Berestetskii modified to ensure positive slope at the origin   *
!   for RG are used.                                                   *
!                                                                      *
!   The arguments are as follows:                                      *
!                                                                      *
!      N    : (Input) The (usual) principal quantum number             *
!      KAPPA: (Input) The relativistic angular quantum number          *
!      Z    : (Input) The effective nuclear charge                     *
!      E    : (Output) The Dirac-Coulomb Eigenenergy                   *
!      RG0  : (Output) Coefficient of the leading term in the          *
!                      series expansion of the large component         *
!                      near the origin                                 *
!      RG   : (Output) r times the large component wavefunction of     *
!                      Akhiezer and Berestetskii                       *
!      RF   : (Output) r times the small component wavefunction of     *
!                      Akhiezer and Berestetskii                       *
!      MTP  : (Output) Maximum tabulation point                        *
!                                                                      *
!   Call(s) to: [LIB92]: CGAMMA.                                       *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last Update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def, only: max_radial_points
      use def_c ,          only: speed_of_light_au, accy
      use grid_c, only: radial_grid, ntp=>radial_point_count
      use tatb_c, only: ta, tb
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use cgamma_i, only: cgamma
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer, intent(in) :: kappa
      integer, intent(out) :: mtp
      real(real64), intent(in) :: z
      real(real64), intent(out) :: e
      real(real64), intent(out) :: rg0
      real(real64), intent(inout) :: rg(max_radial_points)
      real(real64), intent(out) :: rf(max_radial_points)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, nr, nrfac, i, iordr1, iordr2
      real(real64) :: alfa, fn, fkappa, fk, fnr, zalfa, gamma, twogp1, bign, &
         eps, argr, argi, rgamm1, dummy, rgamm2, fac, fg, ff, facn, a, an1, an2&
         , b, bn, fden, bignmk, rho, rhon, f1, f2, ovlfac, cutoff
!-----------------------------------------------
!
!
!
!   Ensure that the principal quantum number is physical
!
      if (n <= 0) then
         write (*, 300)
         write (*, 301) n
         stop
      end if
!
!   Ensure that the angular quantum number is physical
!
      if (kappa == 0) then
         write (*, 300)
         write (*, 302)
         stop
      else if (kappa == n) then
         write (*, 300)
         write (*, 303) kappa, n
         stop
      else if (abs(kappa) > n) then
         write (*, 300)
         write (*, 303) kappa, n
         stop
      end if
!
!   Ensure that the charge is physical
!
      if (z <= 0.0_real64) then
         write (*, 300)
         write (*, 304) z
         stop
      else if (z > speed_of_light_au) then
         write (*, 300)
         write (*, 305) z, speed_of_light_au
         stop
      end if
!
!   Atomic units
!
      alfa = 1.0_real64/speed_of_light_au
!
!   Now determine all the parameters
!
      fn = dble(n)
      fkappa = dble(kappa)
      k = abs(kappa)
      fk = dble(k)
      nr = n - k
      fnr = dble(nr)
      zalfa = z*alfa
      gamma = sqrt(fk*fk - zalfa*zalfa)
      twogp1 = gamma + gamma + 1.0_real64
      bign = sqrt(fn*fn - 2.0_real64*fnr*(fk - gamma))
      eps = 1.0_real64/sqrt(1.0_real64 + (zalfa/(gamma + fnr))**2)
!
!   EPS is the total energy divided by C*C; this must be converted
!   to the units and reference energy of GRASP
!
      e = (1.0_real64 - eps)*speed_of_light_au*speed_of_light_au
!
!   Now the normalization constants
!
      nrfac = 1
      do i = 1, nr
         nrfac = nrfac*i
      end do
!
      argr = twogp1 + fnr
      argi = 0.0_real64
      call cgamma (argr, argi, rgamm1, dummy)
      argr = twogp1
      call cgamma (argr, argi, rgamm2, dummy)
!
      fac = -sqrt(rgamm1)/(rgamm2*sqrt(dble(nrfac)))*sqrt(z/(2.0_real64*bign*bign*(&
         bign - fkappa)))
!
!   Ensure that the slope of the large-component function is
!   positive at the origin
!
      if (kappa > 0) fac = -fac
!
      fg = fac*sqrt(1.0_real64 + eps)
      ff = fac*sqrt(1.0_real64 - eps)
!
!   Now set up the coefficients of the confluent hypergeometric
!   functions  F (-NR+1,2*GAMMA+1;RHO)  and  F (-NR,2*GAMMA+1;RHO)
!   in the workspace arrays  TA  and  TB , respectively
!
      if (nr == 0) then
         iordr1 = 0
         iordr2 = 0
      else
         iordr1 = nr - 1
         iordr2 = nr
      end if
!
      fac = 1.0_real64
      facn = 1.0_real64
      a = -fnr
      an1 = a + 1.0_real64
      an2 = a
      b = twogp1
      bn = b
!
      k = 0
    2 continue
      k = k + 1
      fden = 1.0_real64/(facn*bn)
      if (k <= iordr1) ta(k) = an1*fden
      if (k <= iordr2) then
         tb(k) = an2*fden
         a = a + 1.0_real64
         an1 = an1*(a + 1.0_real64)
         an2 = an2*a
         b = b + 1.0_real64
         bn = bn*b
         fac = fac + 1.0_real64
         facn = facn*fac
         goto 2
      end if
!
!   Now tabulate the function over the entire grid
!
      rg(1) = 0.0_real64
      rf(1) = 0.0_real64
      fac = (z + z)/bign
      bignmk = bign - fkappa
      do i = 2, ntp
         rho = fac*radial_grid(i)
         rhon = rho
         k = 0
         f1 = 1.0_real64
         f2 = 1.0_real64
    3    continue
         k = k + 1
         if (k <= iordr1) f1 = f1 + ta(k)*rhon
         if (k <= iordr2) then
            f2 = f2 + tb(k)*rhon
            rhon = rhon*rho
            goto 3
         end if
         f1 = fnr*f1
         f2 = bignmk*f2
         ovlfac = exp((-0.5_real64*rho))*rho**gamma
         rg(i) = fg*ovlfac*(f1 - f2)
         rf(i) = ff*ovlfac*(f1 + f2)
      end do
!
!   Determine the effective maximum tabulation point based on the
!   cutoff; define the cutoff conservatively
!
      cutoff = accy*0.1_real64
!
      mtp = ntp + 1
    5 continue
      mtp = mtp - 1
!      IF (ABS(RG(MTP)) < CUTOFF) THEN
      if ((abs(rg(mtp)) < cutoff).or.(abs(radial_grid(mtp)) > 1.0e50_real64)) then ! JE: APPLY BOX (R<1.D+50)
         rg(mtp) = 0.0_real64
         rf(mtp) = 0.0_real64
         goto 5
      end if
!
      if (mtp == ntp) write (*, 306) ntp, rg(ntp), cutoff
!
!   Compute the coefficient of R**GAMMA at the origin
!
      rg0 = fg*fac**gamma*(fnr - bignmk)
!
      return
!
  300 format("DCBSRW:")
  301 format(" Principal quantum number is ",1i3)
  302 format(" Angular quantum number is 0")
  303 format(" Angular quantum number (",1i3,") is out of range for",&
         " principal quantum number (",1i3,")")
  304 format(" Nuclear charge (",3p,1d16.7,") is too small")
  305 format(" Nuclear charge (",3p,1d16.7,") exceeds C (",1d16.7,")")
  306 format(/,/,/," ***** Warning in SUBROUTINE DCBSRW *****"/,/,&
         " Radial grid of insufficient extent:"/," P(",1i4,") = ",1p,1d10.3,&
         ", Exceeds cutoff (",1d10.3,")")
      return
!
      end subroutine dcbsrw
end module dcbsrw_i
