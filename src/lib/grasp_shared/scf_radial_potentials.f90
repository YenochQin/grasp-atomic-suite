module scf_radial_potentials
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,    only: max_radial_points
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    ! yp enters the radial Dirac coefficients as yp/r. xp and xq supply the exchange source.
    real(real64) :: yp(max_radial_points), xp(max_radial_points), xq(max_radial_points)
end module scf_radial_potentials
