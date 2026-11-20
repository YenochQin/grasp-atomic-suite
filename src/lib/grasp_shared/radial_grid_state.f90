! Stores the current radial grid and its active point count within the shared workspace capacity.
      module radial_grid_state
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: radial_workspace_capacity
      implicit none (type, external)
      private :: real64
      integer :: radial_point_count
      real(real64) :: radial_grid(radial_workspace_capacity), radial_grid_derivative(radial_workspace_capacity), &
         radial_derivative_over_radius(radial_workspace_capacity)
      real(real64) :: radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      end module radial_grid_state
