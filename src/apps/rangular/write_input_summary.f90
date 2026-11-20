module write_input_summary_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_input_summary
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_input_summary
        ! Module dependencies
        use calculation_file_units, only: ensure_summary_output, summary_output_unit

        !***********************************************************************
        use hamiltonian_symmetry_blocks, only: symmetry_block_count
        use atomic_calculation_state, only: electron_count
        use csf_block_cutoffs, only: active_csf_cutoffs
        use mcp_generation_options, only: diagonal_only, use_first_order_restriction
        use orbital_csf_basis, only: csf_count, orbital_count
        use calendar_clock, only: read_calendar_clock
        use integer_format, only: format_integer
        implicit none (type, external)

        ! Local declarations
        integer :: text_length, symmetry_block_index
        character :: time_text*10, date_text*8
        character :: record*256
        !-----------------------------------------------
        !
        !
        !   Get the date and time of day; make this information the
        !   header of the summary file
        !
        call read_calendar_clock (time_text, date_text)
        call ensure_summary_output()
        write (summary_output_unit, *) "GENMCP run at ", time_text, " on ", date_text, "."
        !
        !   Write out the basic dimensions of the electron cloud
        !
        write (summary_output_unit, *)
        call format_integer (electron_count, record, text_length)
        write (summary_output_unit, *) "There are "//record(1:text_length)//" electrons in the cloud"
        call format_integer (csf_count, record, text_length)
        write (summary_output_unit, *) " in "//record(1:text_length)//" relativistic CSFs"
        call format_integer (orbital_count, record, text_length)
        write (summary_output_unit, *) " based on "//record(1:text_length)//" relativistic subshells."
        !
        !   If the CSFs are not treated uniformly, write out an
        !   informative message
        !
        if (diagonal_only) then
            write (summary_output_unit, *)
            write (summary_output_unit, *) "Only diagonal matrix elements are computed."
        else
            if (use_first_order_restriction) then
                do symmetry_block_index = 1, symmetry_block_count
                    write (summary_output_unit, *)
                    call format_integer (active_csf_cutoffs(symmetry_block_index), record, text_length)
                    write (summary_output_unit, *) "CSFs 1--"//record(1:text_length)//" constitute"//&
                        " the zero-order space."
                end do
            end if
        end if
        !
        return
    end subroutine write_input_summary
end module write_input_summary_interface
