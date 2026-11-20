! Holds energy shifts, eigenvalues, and CSF mixing coefficients for the active, final, and initial states.
module eigenstate_expansions
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: hamiltonian_energy_shift, final_hamiltonian_energy_shift, initial_hamiltonian_energy_shift
    real(real64), pointer :: shifted_eigenvalues(:), final_shifted_eigenvalues(:), initial_shifted_eigenvalues(:)
    real(real64), pointer :: csf_mixing_coefficients(:), final_csf_mixing_coefficients(:), initial_csf_mixing_coefficients(:)
end module eigenstate_expansions
