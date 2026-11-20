module hydrogenic_radial_components
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_radial_points
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: hydrogenic_large_component(max_radial_points), hydrogenic_small_component(max_radial_points)
end module hydrogenic_radial_components
