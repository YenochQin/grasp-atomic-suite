      module grid_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:22:09  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      use parameter_def,   only: radial_workspace_capacity
      implicit none (type, external)
      private :: real64
      integer :: radial_point_count
      real(real64) :: radial_grid(radial_workspace_capacity), radial_grid_derivative(radial_workspace_capacity), &
         radial_derivative_over_radius(radial_workspace_capacity)
      real(real64) :: radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      end module grid_c
