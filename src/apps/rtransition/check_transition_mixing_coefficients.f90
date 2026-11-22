! Check the dimensions and consistency of mixing data for the transition calculation.
module check_transition_mixing_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: check_transition_mixing_coefficients
contains

    !   This routine checks the mixing coefficients
    subroutine check_transition_mixing_coefficients
        ! Module dependencies
        use biorthogonal_orbital_sets, only: final_csf_count, final_orbital_count, initial_csf_count, initial_orbital_count
        use atomic_calculation_state, only: final_electron_count, initial_electron_count
        use eigenstate_expansions, only: final_hamiltonian_energy_shift, initial_hamiltonian_energy_shift, &
            final_shifted_eigenvalues, initial_shifted_eigenvalues, final_csf_mixing_coefficients, initial_csf_mixing_coefficients
        use selected_level_indices, only: final_level_indices, initial_level_indices, final_level_count, initial_level_count, &
            highest_final_level_index, highest_initial_level_index
        use level_symmetry_metadata, only: final_level_parity_signs, initial_level_parity_signs, final_level_two_j_plus_one, &
            initial_level_two_j_plus_one
        implicit none (type, external)

        ! Local declarations
        integer :: i, j

        write (*, *) " ****************"
        write (*, *) " Entering testmix"
        write (*, *) " ****************"
        write (*, *)
        write (*, *) "Initial state"
        write (*, *) "EVALII", (initial_hamiltonian_energy_shift + initial_shifted_eigenvalues(i),i=1,initial_level_count)
        write (*, *) "NELECII,NCFII,NWII,NVECMXII", initial_electron_count, initial_csf_count, initial_orbital_count, &
            highest_initial_level_index
        write (*, *) initial_level_count
        write (*, *) (initial_level_indices(i),i=1,initial_level_count)
        write (*, *) (initial_level_two_j_plus_one(i),initial_level_parity_signs(i),i=1,initial_level_count)
        write (*, *) ((initial_csf_mixing_coefficients(i + (j - 1)*initial_csf_count),i=1,initial_csf_count),j=1, &
            initial_level_count)

        write (*, *) "Final state"
        write (*, *) "EVALFF", (final_hamiltonian_energy_shift + final_shifted_eigenvalues(i),i=1,final_level_count)
        write (*, *) "NELECFF,NCFFF,NWFF,NVECMXFF", final_electron_count, final_csf_count, final_orbital_count, &
            highest_final_level_index
        write (*, *) final_level_count
        write (*, *) (final_level_indices(i),i=1,final_level_count)
        write (*, *) (final_level_two_j_plus_one(i),final_level_parity_signs(i),i=1,final_level_count)
        write (*, *) ((final_csf_mixing_coefficients(i + (j - 1)*final_csf_count),i=1,final_csf_count),j=1,final_level_count)
        write (*, *)
        write (*, *) " ***************"
        write (*, *) " Leaving testmix"
        write (*, *) " ***************"

        return
    end subroutine check_transition_mixing_coefficients
end module check_transition_mixing_coefficients_interface
