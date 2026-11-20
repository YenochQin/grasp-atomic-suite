module nuclear_potential_state
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,    only: max_radial_points
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: first_grid_index_beyond_fermi_c
    ! The stored quantity is -r*V_nuc in atomic units. For a point nucleus it equals Z.
    real(real64) :: negative_radius_times_potential(max_radial_points)
end module nuclear_potential_state
