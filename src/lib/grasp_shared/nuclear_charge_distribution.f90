module nuclear_charge_distribution
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only:  max_radial_points
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: zdist(max_radial_points)
end module nuclear_charge_distribution
