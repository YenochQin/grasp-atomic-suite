module radial_dirac_workspace
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,    only: max_radial_points
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: trial_orbital_last_grid_index
    real(real64) :: trial_large_component(max_radial_points), trial_small_component(max_radial_points)
    real(real64) :: trial_large_origin_parameter, trial_small_origin_parameter
    real(real64) :: large_equation_coefficients(max_radial_points), small_equation_coefficients(max_radial_points), &
        large_radial_source(max_radial_points), small_radial_source(max_radial_points)
end module radial_dirac_workspace
