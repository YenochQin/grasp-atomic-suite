!***********************************************************************
!                                                                      *
      subroutine ncharg
!                                                                      *
!   This routine evaluates the nuclear charge density, and stores it   *
!   in the  common  array  ZDIST .                                     *
!                                                                      *
!   Call(s) to: [LIB92]: ES.                                           *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state,           only: pi,  nuclear_charge, precis
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use nuclear_charge_distribution, only: max_radial_points, zdist
      use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use alternating_exponential_series, only: es
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: radial_index
      real(real64) :: c, a, cba, pi2, abc, abc2, abc3, s2mcba, s3mcba, en, &
         znorm, extrm, zdisti
      logical :: form1, form2
!-----------------------------------------------
!
!
!   Initialize array to zero
!
      zdist(:radial_point_count) = 0.0_real64
!
!   Fermi charge distribution
!
      if (nuclear_model_parameter_count == 2) then
         c = nuclear_model_parameters(1)
         a = nuclear_model_parameters(2)
         cba = c/a
         pi2 = pi*pi
         abc = a/c
         abc2 = abc*abc
         abc3 = abc2*abc
         call es ((-cba), s2mcba, s3mcba)
         en = 1.0_real64 + pi2*abc2 - 6.0_real64*abc3*s3mcba
         znorm = 3.0_real64*nuclear_charge/(4.0_real64*pi*en*c**3)
         form1 = .true.
         form2 = .false.
         do radial_index = 1, radial_point_count
            if (form1) then
               extrm = exp((radial_grid(radial_index)-c)/a)
               zdist(radial_index) = znorm/(1.0_real64 + extrm)
               if (1.0_real64/extrm <= precis) then
                  form1 = .false.
                  form2 = .true.
               end if
            else if (form2) then
               zdisti = znorm*exp((-(radial_grid(radial_index)-c)/a))
               if (abs(zdisti) > 0.0_real64) then
                  zdist(radial_index) = zdisti
               else
                  mtp = radial_index
                  exit
               end if
            end if
         end do
      end if
!
      return
!
      end subroutine ncharg
