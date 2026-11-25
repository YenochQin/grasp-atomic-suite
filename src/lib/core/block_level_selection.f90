! Read the requested ASF indices for each symmetry block and record the largest index needed.
module block_level_selection
    implicit none (type, external)

    ! Module interface
    private
    public :: read_block_level_selection
contains
    subroutine read_block_level_selection(block_labels, abort_handler)
        ! Module dependencies
        use calculation_file_units, only: calculation_log_unit, ensure_calculation_log
        use atomic_calculation_state, only: requested_level_indices, selected_level_count
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, &
            highest_selected_level_by_block, symmetry_block_count
        use standard_io_units, only: stderr_unit, stdin_unit
        use input_error_handling, only: fail_input_io
        use level_selection_parser, only: parse_level_selection, no_levels_selected
        implicit none (type, external)

        ! Arguments
        character(len=8), intent(in) :: block_labels(*)
        interface
            subroutine abort_handler()
                implicit none (type, external)
            end subroutine abort_handler
        end interface
        optional :: abort_handler

        ! Local declarations
        integer :: display_block_index, block_index, block_csf_count, previous_selection_count, input_status, parse_status
        character(len=256) :: selection_input
        character(len=128) :: selection_error

        write (stderr_unit, *) "There are ", symmetry_block_count, " blocks  ", &
            "(block   J/Parity   NCF):"
        write (stderr_unit, "( 4(I3, 1X, A, I8, 5X) )") (display_block_index,trim(block_labels(display_block_index)), &
            & csf_counts_by_block(display_block_index),display_block_index=1&
            ,symmetry_block_count)

        eigenvector_counts_by_block = 0
        highest_selected_level_by_block = 0

        write (stderr_unit, *)
        write (stderr_unit, *) "Enter ASF serial numbers for each block"
        write (stderr_unit, *) "At least one ASF root is required in every J/parity block."

        selected_level_count = 0
        do block_index = 1, symmetry_block_count
            block_csf_count = csf_counts_by_block(block_index)
            write (stderr_unit, *) "Block ", block_index, "   ncf = ", block_csf_count, " id = ", &
                trim(block_labels(block_index))

            previous_selection_count = selected_level_count
            parse_status = no_levels_selected
            read (stdin_unit, "(A)", iostat=input_status) selection_input
            if (input_status == 0) then
                call ensure_calculation_log()
                write (calculation_log_unit, "(A)") trim(selection_input)
                call parse_level_selection (selected_level_count, block_csf_count, selection_input, parse_status)
            end if
            if (parse_status /= 0) then
                write (selection_error, '("Block ",I0," (",A,") requires at least one valid ASF root in [1,",I0,"].")') &
                    block_index, trim(adjustl(block_labels(block_index))), block_csf_count
                if (parse_status == no_levels_selected) selection_error = "Missing ASF selection. "//trim(selection_error)
                call fail_input_io(input_status, stdin_unit, selection_error, abort_handler)
            end if
            eigenvector_counts_by_block(block_index) = selected_level_count - previous_selection_count
            highest_selected_level_by_block(block_index) = &
                maxval(requested_level_indices(previous_selection_count + 1:selected_level_count))
        end do
    end subroutine read_block_level_selection
end module block_level_selection
