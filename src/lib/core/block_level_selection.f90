! Read the requested ASF indices for each symmetry block and record the largest index needed.
module block_level_selection
    implicit none (type, external)

    ! Module interface
    private
    public :: read_block_level_selection
contains
    subroutine read_block_level_selection(block_labels)
        ! Module dependencies
        use calculation_file_units, only: calculation_log_unit, ensure_calculation_log
        use array_memory_management, only: dalloc
        use atomic_calculation_state, only: requested_level_indices, selected_level_count
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, &
            highest_selected_level_by_block, symmetry_block_count
        use standard_io_units, only: stderr_unit
        use level_selection_parser, only: parse_level_selection
        implicit none (type, external)

        ! Arguments
        character (len = 8) :: block_labels(*)

        ! Local declarations
        integer :: display_block_index, selection_index, block_index, block_csf_count, previous_selection_count, &
            & parse_status, highest_selected_level
        character :: selection_input*256

        write (stderr_unit, *) "There are ", symmetry_block_count, " blocks  ", &
            "(block   J/Parity   NCF):"
        write (stderr_unit, "( 4(I3, 1X, A, I8, 5X) )") (display_block_index,trim(block_labels(display_block_index)), &
            & csf_counts_by_block(display_block_index),display_block_index=1&
            ,symmetry_block_count)

        eigenvector_counts_by_block = 0
        highest_selected_level_by_block = 0

        write (stderr_unit, *)
        write (stderr_unit, *) "Enter ASF serial numbers for each block"

        selected_level_count = 0
        do
            do block_index = 1, symmetry_block_count
                block_csf_count = csf_counts_by_block(block_index)
                do
                    write (stderr_unit, *) "Block ", block_index, "   ncf = ", block_csf_count, " id = ", trim(block_labels(&
                        block_index))

                    read (*, "(A)") selection_input
                    call ensure_calculation_log()
                    write(calculation_log_unit,"(a)") trim(selection_input)

                    previous_selection_count = selected_level_count
                    call parse_level_selection (selected_level_count, block_csf_count, selection_input, parse_status)
                    if (selected_level_count == 0) call dalloc (requested_level_indices, "requested_level_indices", &
                        "read_block_level_selection" )
                    if (parse_status >= 0) exit
                end do
                eigenvector_counts_by_block(block_index) = selected_level_count - previous_selection_count

                highest_selected_level = 0
                do selection_index = previous_selection_count + 1, selected_level_count
                    highest_selected_level = max(highest_selected_level,requested_level_indices(selection_index))
                end do
                highest_selected_level_by_block(block_index) = highest_selected_level
            end do

            if (selected_level_count /= 0) exit
            write (stderr_unit, *) "At least one state should be selected"
        end do
    end subroutine read_block_level_selection
end module block_level_selection
