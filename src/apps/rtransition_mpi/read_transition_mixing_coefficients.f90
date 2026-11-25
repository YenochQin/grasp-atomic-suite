module read_transition_mixing_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_transition_mixing_coefficients
contains

    ! Open and read the mixing coefficent files
    !
    ! Written by Per Jonsson
    subroutine read_transition_mixing_coefficients(input_file, init)
        ! Module dependencies
        use mixing_file_format, only: mixing_file, next_mixing_block, read_mixing_levels, read_mixing_vectors
        use array_memory_management, only: alloc
        use biorthogonal_orbital_sets, only: final_csf_count, initial_csf_count
        use eigenstate_expansions, only: final_hamiltonian_energy_shift, initial_hamiltonian_energy_shift, &
            final_shifted_eigenvalues, initial_shifted_eigenvalues, final_csf_mixing_coefficients, initial_csf_mixing_coefficients
        use selected_level_indices, only: final_level_indices, initial_level_indices, final_level_count, initial_level_count
        use level_symmetry_metadata, only: final_level_parity_signs, initial_level_parity_signs, final_level_two_j_plus_one, &
            initial_level_two_j_plus_one
        implicit none (type, external)

        ! Arguments
        type(mixing_file), intent(inout) :: input_file
        integer , intent(in) :: init

        ! Local declarations
        integer :: ib, iatjp, iaspa

        if (init == 1) then
            ! Read the initial state mixing file
            call next_mixing_block(input_file, ib, initial_csf_count, initial_level_count, iatjp, iaspa)
            call alloc (initial_shifted_eigenvalues, initial_level_count, "initial_shifted_eigenvalues", &
                "read_transition_mixing_coefficients")
            call alloc (initial_csf_mixing_coefficients, initial_csf_count*initial_level_count, "initial_csf_mixing_coefficients", &
                "read_transition_mixing_coefficients")
            call alloc (initial_level_indices, initial_level_count,"initial_level_indices", "read_transition_mixing_coefficients" )
            call alloc (initial_level_two_j_plus_one, initial_level_count, "initial_level_two_j_plus_one", &
                "read_transition_mixing_coefficients")
            call alloc (initial_level_parity_signs, initial_level_count, "initial_level_parity_signs", &
                "read_transition_mixing_coefficients")
            initial_level_two_j_plus_one(:initial_level_count) = iatjp
            initial_level_parity_signs(:initial_level_count) = iaspa
            call read_mixing_levels(input_file, initial_level_indices, initial_hamiltonian_energy_shift, &
                initial_shifted_eigenvalues)
            call read_mixing_vectors(input_file, initial_csf_mixing_coefficients, initial_csf_count)

        else
            ! Read the final state mixing file

            call next_mixing_block(input_file, ib, final_csf_count, final_level_count, iatjp, iaspa)
            call alloc (final_shifted_eigenvalues, final_level_count, "final_shifted_eigenvalues", &
                "read_transition_mixing_coefficients")
            call alloc (final_csf_mixing_coefficients, final_csf_count*final_level_count, "final_csf_mixing_coefficients", &
                "read_transition_mixing_coefficients")
            call alloc (final_level_indices, final_level_count,"final_level_indices", "read_transition_mixing_coefficients" )
            call alloc (final_level_two_j_plus_one, final_level_count, "final_level_two_j_plus_one", &
                "read_transition_mixing_coefficients")
            call alloc (final_level_parity_signs, final_level_count, "final_level_parity_signs", &
                "read_transition_mixing_coefficients")
            final_level_two_j_plus_one(:final_level_count) = iatjp
            final_level_parity_signs(:final_level_count) = iaspa
            call read_mixing_levels(input_file, final_level_indices, final_hamiltonian_energy_shift, &
                final_shifted_eigenvalues)
            call read_mixing_vectors(input_file, final_csf_mixing_coefficients, final_csf_count)
            ! Close the initial state mixing  file
        end if

        return
    end subroutine read_transition_mixing_coefficients
end module read_transition_mixing_coefficients_interface
