!***********************************************************************
!                                                                      *
      subroutine radgrd
!                                                                      *
!   This routine sets up the radial grid  R  and the associated arr-   *
!   ays  RP  and  RPOR  in the COMMON block  /GRID/. Different grids   *
!   are generated depending on whether or not  HP  is  zero .          *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use debug_c
      use def_c, only: precis
      use grid_c
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: grid_points_with_padding, radial_index, half_grid_point_count, debug_row_count, debug_row_index, &
         & left_grid_index, right_grid_index
      real(real64) :: exp_grid_step, exp_grid_coordinate, exp_grid_coordinate_minus_one, newton_relative_tolerance, &
         & grid_step_ratio, previous_radius_estimate, radius_estimate, transformed_coordinate, shifted_radius_estimate, &
         & grid_equation_residual, &
         inverse_grid_equation_derivative, radius_correction
!-----------------------------------------------
!
!
!   RPOR(1) is never used in the program: it is arbitrarily
!   set to zero
!
      grid_points_with_padding = radial_point_count + 10
      radial_grid(1) = 0.0_real64
      radial_derivative_over_radius(1) = 0.0_real64
!
!   Now set up the grids
!
      if (asymptotic_grid_step == 0.0_real64) then
!        default comes here
!
!   Exponential grid if HP is zero
!
!   Initializations
!
         radial_grid_derivative(1) = radial_grid_scale
         exp_grid_step = exp(transformed_grid_step)
         exp_grid_coordinate = 1.0_real64
!
!   Set up the arrays R, RP, RPOR
!
         do radial_index = 2, grid_points_with_padding
            exp_grid_coordinate = exp_grid_step*exp_grid_coordinate
            exp_grid_coordinate_minus_one = exp_grid_coordinate - 1.0_real64
            radial_grid(radial_index) = radial_grid_scale*exp_grid_coordinate_minus_one
            radial_grid_derivative(radial_index) = radial_grid_scale*exp_grid_coordinate
            radial_derivative_over_radius(radial_index) = exp_grid_coordinate/exp_grid_coordinate_minus_one
         end do
!
      else
!
!   Asymptotically-linear exponential grid otherwise:
!
!   Initializations
!
         newton_relative_tolerance = 1.0e3_real64*precis
         grid_step_ratio = transformed_grid_step/asymptotic_grid_step
         radial_grid_derivative(1) = radial_grid_scale/(grid_step_ratio*radial_grid_scale + 1.0_real64)
         previous_radius_estimate = 0.0_real64
         radius_estimate = 0.0_real64
!
!   Set up the arrays R, RP, RPOR
!
         do radial_index = 2, grid_points_with_padding
!
            transformed_coordinate = transformed_grid_step*dble(radial_index - 1)
!
!   Solve the implicit equation for R using the Newton-Raphson
!   method
!
    2       continue
            shifted_radius_estimate = radius_estimate + radial_grid_scale
            grid_equation_residual = log(shifted_radius_estimate/radial_grid_scale) + grid_step_ratio*radius_estimate &
               & - transformed_coordinate
            inverse_grid_equation_derivative = shifted_radius_estimate/(grid_step_ratio*shifted_radius_estimate + 1.0_real64)
            radius_correction = -grid_equation_residual*inverse_grid_equation_derivative
            radius_estimate = previous_radius_estimate + radius_correction
!
            if (abs(radius_correction/radius_estimate) < newton_relative_tolerance) then
               radial_grid(radial_index) = radius_estimate
               shifted_radius_estimate = radius_estimate + radial_grid_scale
               inverse_grid_equation_derivative = shifted_radius_estimate/(grid_step_ratio*shifted_radius_estimate + 1.0_real64)
               radial_grid_derivative(radial_index) = inverse_grid_equation_derivative
               radial_derivative_over_radius(radial_index) = inverse_grid_equation_derivative/radius_estimate
            else
               previous_radius_estimate = radius_estimate
               go to 2
            endif
!
         end do
!
      endif
!
!   Debug printout
!
      if (ldbpr(1)) then
         write (99, 300)
         half_grid_point_count = radial_point_count/2
         if (2*half_grid_point_count == radial_point_count) then
            debug_row_count = half_grid_point_count
         else
            debug_row_count = half_grid_point_count + 1
         endif
         do debug_row_index = 1, debug_row_count
            left_grid_index = debug_row_index
            right_grid_index = left_grid_index + debug_row_count
            if (right_grid_index <= radial_point_count) then
               write (99, 301) radial_grid(left_grid_index), radial_grid_derivative(left_grid_index), &
                  & radial_derivative_over_radius(left_grid_index), radial_grid(right_grid_index), &
                  & radial_grid_derivative(right_grid_index), &
                  radial_derivative_over_radius(right_grid_index)
            else if (left_grid_index <= radial_point_count) then
               write (99, 301) radial_grid(left_grid_index), radial_grid_derivative(left_grid_index), &
                  & radial_derivative_over_radius(left_grid_index)
            endif
         end do
      endif
!
      return
!
  300 format(/,'From SUBROUTINE RADGRD:'/,2(&
         ' -------- r -------- -------- r'' -------',' ------- r''/r ------'))
  301 format(1p,6(1x,1d19.12))
      return
!
      end subroutine radgrd
