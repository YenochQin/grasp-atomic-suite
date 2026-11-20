! Stores the current radial grid and its active point count within the shared workspace capacity.
module radial_grid_state
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: radial_workspace_capacity
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: radial_point_count
    ! Points are uniform in t. The arrays store r(t), dr/dt, and (dr/dt)/r with the last set to zero at r=0.
    real(real64) :: radial_grid(radial_workspace_capacity), radial_grid_derivative(radial_workspace_capacity), &
        radial_derivative_over_radius(radial_workspace_capacity)
    real(real64) :: radial_grid_scale, transformed_grid_step, asymptotic_grid_step
end module radial_grid_state
