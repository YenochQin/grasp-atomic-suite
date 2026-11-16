      module grid_c
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  06:22:09  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      use parameter_def,   only: radial_workspace_capacity
      integer :: radial_point_count
      real(double), dimension(radial_workspace_capacity) :: radial_grid, radial_grid_derivative, radial_derivative_over_radius
      real(double) :: radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      end module grid_c
