! Stores orbital quantum numbers, energies, and occupations for the CSF basis.
module orbital_csf_basis
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: iso_fortran_env, only: int8
    use grasp_parameters,   only:  max_orbitals
    implicit none (type, external)

    ! Module interface
    private :: real64, int8

    ! Module state
    character(len=2) :: orbital_symmetry_labels(max_orbitals)
    character(len=2) :: left_orbital_labels(max_orbitals)
    ! The radial solvers use a positive binding-energy convention for bound orbitals.
    real(real64) :: orbital_energy_parameters(max_orbitals), orbital_origin_exponents(max_orbitals), &
        previous_orbital_energy_changes(max_orbitals)
    integer :: csf_count, orbital_count, left_csf_count, left_orbital_count
    ! Occupation storage is indexed by orbital, then by CSF.
    integer(int8), allocatable :: csf_subshell_occupations(:, :)
    ! The paired left basis uses default-kind occupations in the same index order.
    integer, allocatable :: left_occupations(:, :)
    integer :: orbital_principal_n(max_orbitals), orbital_kappa(max_orbitals), left_principal_n(max_orbitals), &
        left_kappa(max_orbitals)
    integer :: orbital_l(max_orbitals), orbital_two_j(max_orbitals)
end module orbital_csf_basis
