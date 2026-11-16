module start_i
    implicit none (type, external)
    private
    public :: start
contains
!***********************************************************************
!                                                                      *
      subroutine start(orbital_index, radial_equation_type, origin_slope_parameter, large_component_values, &
         & small_component_origin_coefficient, small_component_values)
!                                                                      *
!   This subroutine sets up  P(1:6), Q(1:6),  required  to start the   *
!   integration for programs  OUT  and  SBSTEP .                       *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!      IORB : (Input) Index of the orbital                             *
!      ITYPE: (Input) 1 = homogeneous equation; 2 = inhomogeneous      *
!             equation; 3 = variational equation                       *
!      P0   : (Input) Slope parameter                                  *
!      P    : (Output) P(1:6) are tabulated by this program            *
!      Q0   : (Output) First term in the series expansion of Q         *
!      Q    : (Output) Q(1:6) are tabulated by this program            *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_radial_points
      use cnc_c, only: cnc5c, cnc6c, double
      use def_c,           only: speed_of_light_au, nuclear_charge, accy
      use grid_c
      use npar_c, only: double, nuclear_model_parameter_count, nuclear_model_parameters
      use orb_c
      use pote_c,          only: yp, xp, xq
      use scf_c,           only: ndcof, nda, da
      use wave_c,          only: orbital_origin_parameters,radial_large_components,radial_small_components
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: orbital_index
      integer, intent(in) :: radial_equation_type
      real(real64), intent(in) :: origin_slope_parameter
      real(real64), intent(out) :: small_component_origin_coefficient
      real(real64), intent(inout) :: large_component_values(max_radial_points)
      real(real64), intent(inout) :: small_component_values(max_radial_points)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, parameter :: mxiter = 36
      integer :: kiorb, i, jorb, niter, j
      real(real64) :: rdp(6), rdq(6), rsep(6), rseq(6)
      real(real64) :: spest(2:6), sqest(2:6)
      real(real64) :: obc, zonc, giorb, fkiorb, omgi, omgmk, omgpk, rsep1, &
         rseq1, pzero, p1, q1, sump, sumq, factor, csq, twocsq, energy, enefac&
         , ri, rpi, rirpi, ypirpi, difmaw, difmax, coefij, pi, qi, rj
!-----------------------------------------------
!
!
!   Initialization
!
      obc = 1.0_real64/speed_of_light_au
      zonc = nuclear_charge*obc
      giorb = orbital_origin_exponents(orbital_index)
      kiorb = orbital_kappa(orbital_index)
      fkiorb = dble(kiorb)
      omgi = 1.0_real64 - giorb
!
      omgmk = omgi - fkiorb
      omgpk = omgi + fkiorb
!
!   Determine P(1), Q(1): THESE STORE  R**(-GAMMA)*(P(1),Q(1));
!   set up  RSEP  and  RSEQ , the inhomogeneous terms
!
      if (radial_equation_type==1 .or. radial_equation_type==2) then
         if (nuclear_model_parameter_count == 0) then
            large_component_values(1) = origin_slope_parameter
            if (kiorb < 0) then
               small_component_values(1) = -origin_slope_parameter*zonc/(giorb - fkiorb)
            else
               small_component_values(1) = origin_slope_parameter*(giorb + fkiorb)/zonc
            end if
         else
            if (kiorb < 0) then
               large_component_values(1) = origin_slope_parameter
               small_component_values(1) = 0.0_real64
            else
               large_component_values(1) = 0.0_real64
               small_component_values(1) = origin_slope_parameter*(giorb + fkiorb)/zonc
            end if
         end if
         if (radial_equation_type == 1) then
            rsep = 0.0_real64
            rseq = 0.0_real64
         else
            rsep1 = 0.0_real64
            rseq1 = 0.0_real64
            do i = 1, ndcof
               jorb = nda(i)
               pzero = orbital_origin_parameters(jorb)
               if (nuclear_model_parameter_count == 0) then
                  p1 = pzero
                  if (kiorb < 0) then
                     q1 = -pzero*zonc/(giorb - fkiorb)
                  else
                     q1 = pzero*(giorb + fkiorb)/zonc
                  end if
                  sump = zonc*q1
                  sumq = -zonc*p1
               else
                  if (kiorb < 0) then
                     p1 = pzero
                     q1 = 0.0_real64
                  else
                     p1 = 0.0_real64
                     q1 = pzero*(giorb + fkiorb)/zonc
                  end if
                  sump = 0.0_real64
                  sumq = 0.0_real64
               end if
               factor = da(i)
               rsep1 = rsep1 + factor*(sump + omgmk*p1)
               rseq1 = rseq1 + factor*(sumq + omgpk*q1)
            end do
            factor = radial_grid_derivative(1)
            rsep(1) = factor*rsep1
            rseq(1) = factor*rseq1
            do i = 2, 6
               factor = -radial_grid_derivative(i)*radial_grid(i)**(-giorb)
               rsep(i) = factor*xp(i)
               rseq(i) = factor*xq(i)
            end do
         end if
      else if (radial_equation_type == 3) then
         large_component_values(1) = 0.0_real64
         small_component_values(1) = 0.0_real64
         rsep(1) = 0.0_real64
         rseq(1) = 0.0_real64
         do i = 2, 6
            factor = obc*radial_grid_derivative(i)*radial_grid(i)**omgi
            rsep(i) = -factor*radial_small_components(i,orbital_index)
            rseq(i) = factor*radial_large_components(i,orbital_index)
         end do
      end if
      small_component_origin_coefficient = small_component_values(1)
!
!   Set up  RDP  and  RDQ
!
      csq = speed_of_light_au*speed_of_light_au
      twocsq = csq + csq
      energy = orbital_energy_parameters(orbital_index)
      enefac = twocsq - energy
      do i = 1, 6
         ri = radial_grid(i)
         rpi = radial_grid_derivative(i)
         rirpi = ri*rpi
         ypirpi = yp(i)*rpi
         rdp(i) = -obc*(enefac*rirpi + ypirpi)
         rdq(i) = -obc*(energy*rirpi - ypirpi)
      end do
!
!   Determine  P(2:6) , Q(2:6)
!
!   Initilizations for the iterations
!
      niter = 0
      p1 = large_component_values(1)
      q1 = small_component_values(1)
      difmaw = max(abs(p1),abs(q1))
!
      large_component_values(2:6) = p1
      small_component_values(2:6) = q1
!
!   This is the largest factor by which any result will be
!   multiplied
!
      factor = radial_grid(6)**giorb
!
!   Now iterate
!
    7 continue
      niter = niter + 1
      difmax = 0.0_real64
      do j = 2, 6
         sump = sum(cnc6c(:,j)*(omgmk*radial_grid_derivative(:6)*large_component_values(:6)-rdp*small_component_values(:6)+rsep))
         sumq = sum(cnc6c(:,j)*(omgpk*radial_grid_derivative(:6)*small_component_values(:6)-rdq*large_component_values(:6)+rseq))
         rj = radial_grid(j)
         sump = sump/rj
         sumq = sumq/rj
         spest(j) = sump
         sqest(j) = sumq
         difmax = max(difmax,abs(sump - large_component_values(j)))
         difmax = max(difmax,abs(sumq - small_component_values(j)))
      end do
!zou  IF (DIFMAX .LT. DIFMAW) THEN
      large_component_values(2:6) = spest(:6)
      small_component_values(2:6) = sqest(:6)
      difmaw = difmax
      difmax = difmax*factor
      if (difmax > accy) then
         if (niter < mxiter) then
            goto 7
         end if
      end if
!     ELSE
!        DIFMAX = DIFMAX*FACTOR
!        IF (DIFMAX .GT. ACCY) THEN
!           WRITE (*,300) NP(IORB),NH(IORB),DIFMAX,NITER,ACCY
!        ENDIF
!zou  ENDIF
! not convergent, using the initial P,Q
      if (difmax > accy) then
         large_component_values(2:6) = p1
         small_component_values(2:6) = q1
      end if
!zou
!
!   All done
!
!   This is always true in GRASP2
!
      large_component_values(1) = 0.0_real64
      small_component_values(1) = 0.0_real64
!
      do i = 2, 6
         factor = radial_grid(i)**giorb
         large_component_values(i) = factor*large_component_values(i)
         small_component_values(i) = factor*small_component_values(i)
      end do
!
      return
!
  300 format("START: ",1i2,1a2," subshell: accuracy ",1p,1d8.1,/,&
         " attained after ",1i2," iterations; this fails the"/,&
         " accuracy criterion ",d8.1,".")
      return
!
      end subroutine start
end module start_i
