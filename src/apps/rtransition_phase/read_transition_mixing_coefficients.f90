module read_transition_mixing_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_transition_mixing_coefficients
contains

    !***********************************************************************
    !                                                                      *
    subroutine read_transition_mixing_coefficients(name, inpci, init)
        ! Module dependencies
        use transition_matrix_workspace, only: final_mixing_unit, initial_mixing_unit
        !                                                                      *
        !   Open and read the mixing coefficent files                          *
        !                                                                      *
        !   Written by Per Jonsson                                             *
        !                                                                      *
        !***********************************************************************
        use array_memory_management, only: alloc
        use biorthogonal_orbital_sets, only: final_csf_count, initial_csf_count
        use eigenstate_expansions, only: final_hamiltonian_energy_shift, initial_hamiltonian_energy_shift, &
            final_shifted_eigenvalues, initial_shifted_eigenvalues, final_csf_mixing_coefficients, initial_csf_mixing_coefficients
        use selected_level_indices, only: final_level_indices, initial_level_indices, final_level_count, initial_level_count
        use level_symmetry_metadata, only: final_level_parity_signs, initial_level_parity_signs, final_level_two_j_plus_one, &
            initial_level_two_j_plus_one
        implicit none (type, external)

        ! Arguments
        integer  :: inpci
        integer , intent(in) :: init
        character  :: name(2)*24

        ! Local declarations
        integer :: ib, iatjp, iaspa, i, j
        !-----------------------------------------------
        !
        if (init == 1) then
            !
            !   Read the initial state mixing file
            !
            read (initial_mixing_unit) ib, initial_csf_count, initial_level_count, iatjp, iaspa
            call alloc (initial_shifted_eigenvalues, initial_level_count, "initial_shifted_eigenvalues", &
                "read_transition_mixing_coefficients")
            call alloc (initial_csf_mixing_coefficients, initial_csf_count*initial_level_count, "initial_csf_mixing_coefficients", &
                "read_transition_mixing_coefficients")
            call alloc (initial_level_indices, initial_level_count,"initial_level_indices", "read_transition_mixing_coefficients" )
            call alloc (initial_level_two_j_plus_one, initial_level_count, "initial_level_two_j_plus_one", &
                "read_transition_mixing_coefficients")
            call alloc (initial_level_parity_signs, initial_level_count, "initial_level_parity_signs", &
                "read_transition_mixing_coefficients")
            read (initial_mixing_unit) (initial_level_indices(i),i=1,initial_level_count)
            initial_level_two_j_plus_one(:initial_level_count) = iatjp
            initial_level_parity_signs(:initial_level_count) = iaspa
            read (initial_mixing_unit) initial_hamiltonian_energy_shift, (initial_shifted_eigenvalues(i),i=1,initial_level_count)

            read (initial_mixing_unit) ((initial_csf_mixing_coefficients(i + (j - 1)*initial_csf_count),i=1, &
                initial_csf_count),j=1,initial_level_count)
        !
        !     CLOSE(68)
        !
        else
            !
            !   Read the final state mixing file
            !

            read (final_mixing_unit) ib, final_csf_count, final_level_count, iatjp, iaspa
            call alloc (final_shifted_eigenvalues, final_level_count, "final_shifted_eigenvalues", &
                "read_transition_mixing_coefficients")
            call alloc (final_csf_mixing_coefficients, final_csf_count*final_level_count, "final_csf_mixing_coefficients", &
                "read_transition_mixing_coefficients")
            call alloc (final_level_indices, final_level_count,"final_level_indices", "read_transition_mixing_coefficients" )
            call alloc (final_level_two_j_plus_one, final_level_count, "final_level_two_j_plus_one", &
                "read_transition_mixing_coefficients")
            call alloc (final_level_parity_signs, final_level_count, "final_level_parity_signs", &
                "read_transition_mixing_coefficients")
            read (final_mixing_unit) (final_level_indices(i),i=1,final_level_count)
            final_level_two_j_plus_one(:final_level_count) = iatjp
            final_level_parity_signs(:final_level_count) = iaspa
            read (final_mixing_unit) final_hamiltonian_energy_shift, (final_shifted_eigenvalues(i),i=1,final_level_count)

            read (final_mixing_unit) ((final_csf_mixing_coefficients(i + (j - 1)*final_csf_count),i=1,final_csf_count),j=1, &
                final_level_count)
            !
            !   Close the initial state mixing  file
            !
            !      CLOSE(78)
        end if

        return
    end subroutine read_transition_mixing_coefficients
end module read_transition_mixing_coefficients_interface
