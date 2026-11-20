module bessel_function_workspace
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only:  max_radial_points
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: breit_bessel_frequencies(2)
    real(real64) :: breit_spherical_bessel_values(2, 2, max_radial_points), breit_spherical_neumann_values(2, 2, max_radial_points)
    real(real64) :: transition_spherical_bessel_values(max_radial_points, 3)
    real(real64) :: cross_component_product_difference(max_radial_points), same_component_product_sum(max_radial_points)
end module bessel_function_workspace
