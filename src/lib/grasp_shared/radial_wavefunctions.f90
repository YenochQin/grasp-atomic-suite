module radial_wavefunctions
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only:  max_orbitals
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    ! The last usable sample can differ by orbital. The origin coefficient is stored separately
    ! from the sampled radial components, which vanish at the origin.
    integer :: orbital_last_grid_index(max_orbitals)
    real(real64) :: orbital_origin_parameters(max_orbitals)
    ! Retain separate initial and final orbitals for transition calculations.
    integer :: final_orbital_last_grid_index(max_orbitals)
    real(real64) :: final_orbital_origin_parameters(max_orbitals)
    integer :: initial_orbital_last_grid_index(max_orbitals)
    real(real64) :: initial_orbital_origin_parameters(max_orbitals)
    ! Radial components are indexed by grid point, then by orbital.
    real(real64), pointer :: radial_large_components(:, :),radial_small_components(:, :)
    real(real64), pointer :: final_radial_large_components(:, :),final_radial_small_components(:, :), &
        initial_radial_large_components(:, :),initial_radial_small_components(:, :)
end module radial_wavefunctions
