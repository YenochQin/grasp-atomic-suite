module scf_radial_potentials
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,    only: max_radial_points
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    ! Divide the orbital Y potential by radius in the radial Dirac equations.
    ! The exchange source is stored separately for the large and small components.
    real(real64) :: orbital_y_potential(max_radial_points), large_exchange_source(max_radial_points), &
        small_exchange_source(max_radial_points)
end module scf_radial_potentials
