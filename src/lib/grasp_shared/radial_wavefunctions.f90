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
    ! The ff and ii arrays retain separate final and initial orbitals for transitions.
    integer :: mfff(max_orbitals)
    real(real64) :: pzff(max_orbitals)
    integer :: mfii(max_orbitals)
    real(real64) :: pzii(max_orbitals)
    ! Radial components are indexed by grid point, then by orbital.
    real(real64), pointer :: radial_large_components(:, :),radial_small_components(:, :)
    real(real64), pointer :: pfff(:, :),qfff(:, :),pfii(:, :),qfii(:, :)
end module radial_wavefunctions
