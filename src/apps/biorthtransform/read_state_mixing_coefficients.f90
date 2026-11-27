module read_state_mixing_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_state_mixing_coefficients
contains

    ! Open the state mixing stream and read the next symmetry block.
    !
    ! Written by Farid A. Parpia. Last revision: 25 Dec 1992
    subroutine read_state_mixing_coefficients(name, inpci, iblk, input_file)
        ! Module dependencies
        use csf_descriptor_reader, only: loaded_csf_basis_id
        use mixing_file_format, only: mixing_file, open_mixing_input, next_mixing_block, read_mixing_levels, read_mixing_vectors
        use calculation_state_files, only: mixing_input_path
        use, intrinsic :: iso_fortran_env, only: output_unit
        use array_memory_management, only: alloc
        use atomic_calculation_state, only: electron_count
        use orbital_csf_basis, only: csf_count, orbital_count
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
        use selected_level_indices, only: active_level_indices, active_level_count
        use level_symmetry_metadata, only: level_parity_signs, level_two_j_plus_one
        use transition_basis_blocks, only: transition_block_count
        use input_error_handling, only: fail_input_io
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: inpci
        integer, intent(in) :: iblk
        character, intent(in) :: name*24

        type(mixing_file), intent(inout) :: input_file

        ! Local declarations
        integer :: nb, iatjp, iaspa

        if (iblk == 1) then
            call open_mixing_input(input_file, mixing_input_path(name, inpci == 0), loaded_csf_basis_id)
            if (electron_count /= input_file%electron_count .or. orbital_count /= input_file%orbital_count .or. &
                transition_block_count /= input_file%block_count) &
                call fail_input_io(0, input_file%unit, "Mixing stream dimensions do not match the biorthogonal basis")
        end if
        ! Load data from the  rscf.mix  file
        write (output_unit, *) "Loading MIXing Coefficients File ..."
        call next_mixing_block(input_file, nb, csf_count, active_level_count, iatjp, iaspa)
        if (nb /= iblk .or. active_level_count == 0) error stop "Cannot biorthogonally transform an empty or unordered block"
        allocate(shifted_eigenvalues(active_level_count))
        allocate(csf_mixing_coefficients(csf_count*active_level_count))
        call alloc (active_level_indices, active_level_count, "active_level_indices", "read_state_mixing_coefficients")
        call alloc (level_two_j_plus_one, active_level_count, "level_two_j_plus_one", "read_state_mixing_coefficients")
        call alloc (level_parity_signs, active_level_count, "level_parity_signs", "read_state_mixing_coefficients")
        ! These arrays are deallocated in mcp
        level_two_j_plus_one(:active_level_count) = iatjp
        level_parity_signs(:active_level_count) = iaspa
        call read_mixing_levels(input_file, active_level_indices, hamiltonian_energy_shift, shifted_eigenvalues)
        call read_mixing_vectors(input_file, csf_mixing_coefficients, csf_count)
        write (output_unit, *) " ... load complete;"
        return
    end subroutine read_state_mixing_coefficients
end module read_state_mixing_coefficients_interface
