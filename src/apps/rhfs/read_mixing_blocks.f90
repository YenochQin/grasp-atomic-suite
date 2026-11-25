! Load level energies and mixing coefficients from the selected block-structured mixing file.
module read_mixing_blocks_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_mixing_blocks
contains

    !    Reads mixing coefficient file from block-structured format
    !
    ! Note:
    !     eav is not compatible with the non-block version if some blocks
    !     were not diagonalized
    !
    !     This is a modified version of cvtmix.f
    ! written by Per Jonsson, September 2003
    subroutine read_mixing_blocks(state_name, ci_input_mode)
        ! Module dependencies
        use mixing_file_format, only: mixing_file, open_mixing_input, next_mixing_block, read_mixing_levels, read_mixing_vectors
        use csf_descriptor_reader, only: loaded_csf_basis_id
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use calculation_state_files, only: mixing_input_path
        use array_memory_management, only: alloc
        use atomic_calculation_state, only: electron_count
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
        use orbital_csf_basis, only: csf_count, orbital_count
        use selected_level_indices, only: active_level_indices, active_level_count
        use level_symmetry_metadata, only: level_parity_signs, level_two_j_plus_one
        implicit none (type, external)

        ! Arguments and local declarations
        integer, intent(in)           :: ci_input_mode
        character(len=24), intent(in) :: state_name
        integer :: total_csf_count, total_level_count, block_count, &
                & preceding_level_count, &
            preceding_csf_count, preceding_mixing_size, contributing_csf_count, block_index, stored_block_index, &
                & block_csf_count, block_level_count, block_two_j_plus_one, block_parity_sign
        real(real64) :: weighted_energy_shift_sum
        character :: file_path*256
        type(mixing_file) :: input_file

        file_path = mixing_input_path(state_name, ci_input_mode == 0)
        call open_mixing_input(input_file, file_path, loaded_csf_basis_id)
        electron_count = input_file%electron_count
        orbital_count = input_file%orbital_count
        total_csf_count = input_file%csf_count
        total_level_count = input_file%level_count
        block_count = input_file%block_count
        write (*, *) "   nelec  = ", electron_count
        write (*, *) "   ncftot = ", total_csf_count
        write (*, *) "   nw     = ", orbital_count
        write (*, *) "   nblock = ", block_count
        write (*, *)

        ! Allocate memory for old format data

        call alloc (shifted_eigenvalues, total_level_count, "shifted_eigenvalues", "read_mixing_blocks")
        call alloc (csf_mixing_coefficients, total_csf_count*total_level_count, "csf_mixing_coefficients", "read_mixing_blocks")
        call alloc (active_level_indices, total_level_count, "active_level_indices", "read_mixing_blocks")
        call alloc (level_two_j_plus_one, total_level_count, "level_two_j_plus_one", "read_mixing_blocks")
        call alloc (level_parity_signs, total_level_count, "level_parity_signs", "read_mixing_blocks")

        ! Initialize mixing coefficients to zero; others are fine
        csf_mixing_coefficients(:total_level_count*total_csf_count) = 0.0_real64

        ! Initialize counters and sum registers
        !
        !    nvecpat:    total number of eigenstates of the previous blocks
        !    ncfpat:     total number of CSF of the previous blocks
        !    nvecsizpat: vector size of the previous blocks
        !    eavsum:     sum of diagonal elements of the previous blocks where
        !                at least one eigenstate is calculated
        !    neavsum:    total number CSF contributing to eavsum

        preceding_level_count = 0
        preceding_csf_count = 0
        preceding_mixing_size = 0
        contributing_csf_count = 0
        weighted_energy_shift_sum = 0.0_real64

        write (*, *) "  block     ncf     nev    2j+1  parity"
        do block_index = 1, block_count

            call next_mixing_block(input_file, stored_block_index, block_csf_count, block_level_count, &
                block_two_j_plus_one, block_parity_sign)
            write (*, "(5I8)") stored_block_index, block_csf_count, block_level_count, block_two_j_plus_one, block_parity_sign
            if (block_index /= stored_block_index) stop "jb .NE. nb"

            if (block_level_count > 0) then

                call read_mixing_levels(input_file, &
                    active_level_indices(preceding_level_count+1:preceding_level_count+block_level_count), &
                    hamiltonian_energy_shift, shifted_eigenvalues(preceding_level_count+1:preceding_level_count+block_level_count))
                ! ivec(i)   = ivec(i) + ncfpat ! serial # of the state
                level_two_j_plus_one(preceding_level_count+1:block_level_count+preceding_level_count) = block_two_j_plus_one
                level_parity_signs(preceding_level_count+1:block_level_count+preceding_level_count) = block_parity_sign

                ! ...Construct the true energy by adding up the average
                shifted_eigenvalues(preceding_level_count+1:block_level_count+ &
                & preceding_level_count) = shifted_eigenvalues(preceding_level_count+1:block_level_count+preceding_level_count) + &
                    hamiltonian_energy_shift
                ! ...For overal (all blocks) average energy
                weighted_energy_shift_sum = weighted_energy_shift_sum + hamiltonian_energy_shift*block_csf_count
                contributing_csf_count = contributing_csf_count + block_csf_count

                call read_mixing_vectors(input_file, &
                    csf_mixing_coefficients(preceding_mixing_size+preceding_csf_count+1:), total_csf_count)
            end if

            preceding_level_count = preceding_level_count + block_level_count
            preceding_csf_count = preceding_csf_count + block_csf_count
            preceding_mixing_size = preceding_mixing_size + block_level_count*total_csf_count

        end do

        ! ...Here eav is the average energy of the blocks where at least
        !    one eigenstate is calculated. It is not the averge of the
        !    total Hamiltonian.

        hamiltonian_energy_shift = weighted_energy_shift_sum/contributing_csf_count

        if (total_csf_count /= contributing_csf_count) write (output_unit, *) &
            "Not all blocks are diagonalized --- Average E ", "not correct"

        ! ...Substrct the overal average energy
        shifted_eigenvalues(:total_level_count) = shifted_eigenvalues(:total_level_count) - hamiltonian_energy_shift

        close(input_file%unit)

        csf_count = total_csf_count
        active_level_count = total_level_count

        return
    end subroutine read_mixing_blocks
end module read_mixing_blocks_interface
