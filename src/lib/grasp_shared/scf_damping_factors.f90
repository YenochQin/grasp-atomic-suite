! Stores separate SCF damping factors for radial orbitals and CSF mixing coefficients.
module scf_damping_factors
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,    only: max_orbitals
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: orbital_damping_factors(max_orbitals)
    real(real64), pointer :: mixing_damping_factors(:)
end module scf_damping_factors
