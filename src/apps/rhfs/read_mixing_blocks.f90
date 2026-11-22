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
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use calculation_state_files, only: mixing_input_path
        use array_memory_management, only: alloc
        use atomic_calculation_state, only: electron_count
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
        use orbital_csf_basis, only: csf_count, orbital_count
        use selected_level_indices, only: active_level_indices, active_level_count
        use level_symmetry_metadata, only: level_parity_signs, level_two_j_plus_one
        use standard_io_units, only: stderr_unit
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: mixing_input_unit
        integer, intent(in)           :: ci_input_mode
        character(len=24), intent(in) :: state_name
        integer :: open_status, io_status, total_csf_count, total_level_count, vector_size, block_count, item_index, &
                & preceding_level_count, &
            preceding_csf_count, preceding_mixing_size, contributing_csf_count, block_index, stored_block_index, &
                & block_csf_count, block_level_count, block_two_j_plus_one, block_parity_sign, block_level_index
        real(real64) :: weighted_energy_shift_sum
        character :: file_path*256, form*11, mixing_file_signature*6, status*3
        ! The  .mix  file is UNFORMATTED; it must exist
        file_path = mixing_input_path(state_name, ci_input_mode == 0)
        form = "UNFORMATTED"
        status = "OLD"
        call open_calculation_file (mixing_input_unit, file_path, form, status, open_status)
        if (open_status == 1) then
            write (stderr_unit, *) "Error when opening", file_path
            stop
        end if
        ! Check the header of the file; if not as expected, try again
        read (mixing_input_unit, iostat=io_status) mixing_file_signature
        if (io_status/=0 .or. mixing_file_signature/="G92MIX") then
            write (stderr_unit, *) "Not a GRASP92 MIXing Coefficients File;"
            close(mixing_input_unit)
            stop
        end if

        read (mixing_input_unit) electron_count, total_csf_count, orbital_count, total_level_count, vector_size, block_count
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

            read (mixing_input_unit) stored_block_index, block_csf_count, block_level_count, block_two_j_plus_one, block_parity_sign
            write (*, "(5I8)") stored_block_index, block_csf_count, block_level_count, block_two_j_plus_one, block_parity_sign
            if (block_index /= stored_block_index) stop "jb .NE. nb"

            if (block_level_count > 0) then

                read (mixing_input_unit) (active_level_indices(preceding_level_count + item_index),item_index=1,block_level_count)
                ! ivec(i)   = ivec(i) + ncfpat ! serial # of the state
                level_two_j_plus_one(preceding_level_count+1:block_level_count+preceding_level_count) = block_two_j_plus_one
                level_parity_signs(preceding_level_count+1:block_level_count+preceding_level_count) = block_parity_sign

                read (mixing_input_unit) hamiltonian_energy_shift, (shifted_eigenvalues(preceding_level_count+item_index), &
                & item_index=1,block_level_count)

                ! ...Construct the true energy by adding up the average
                shifted_eigenvalues(preceding_level_count+1:block_level_count+ &
                & preceding_level_count) = shifted_eigenvalues(preceding_level_count+1:block_level_count+preceding_level_count) + &
                    hamiltonian_energy_shift
                ! ...For overal (all blocks) average energy
                weighted_energy_shift_sum = weighted_energy_shift_sum + hamiltonian_energy_shift*block_csf_count
                contributing_csf_count = contributing_csf_count + block_csf_count

                read (mixing_input_unit) ((csf_mixing_coefficients(preceding_mixing_size+preceding_csf_count+item_index+ &
                & (block_level_index-1)*total_csf_count),item_index=1,block_csf_count),block_level_index=1,&
                    block_level_count)
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

        close(mixing_input_unit)

        csf_count = total_csf_count
        active_level_count = total_level_count

        return
    end subroutine read_mixing_blocks
end module read_mixing_blocks_interface
