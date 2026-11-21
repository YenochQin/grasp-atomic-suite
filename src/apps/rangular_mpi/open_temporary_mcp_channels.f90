module open_temporary_mcp_channels_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_temporary_mcp_channels
contains

    !***********************************************************************
    !                                                                      *
    subroutine open_temporary_mcp_channels(symmetry_block_index, max_coulomb_multipole_rank, output_file_prefix)
        ! Module dependencies
        use mcp_generation_options, only: mcp_temporary_units
        !***********************************************************************
        use integer_format, only: format_integer
        use file_open, only: open_calculation_file
        use standard_io_units, only: stderr_unit
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: symmetry_block_index
        integer , intent(in) :: max_coulomb_multipole_rank
        character , intent(in) :: output_file_prefix*(*)

        ! Local declarations
        integer :: mcp_channel_index, text_length, open_status, opened_channel_index, prefix_length
        character :: number_text*2
        !-----------------------------------------------
        !
        ! All files  filehead.XX  are UNFORMATTED;
        !
        prefix_length = len_trim(output_file_prefix)
        do mcp_channel_index = 30, 32 + max_coulomb_multipole_rank
            call format_integer (mcp_channel_index, number_text, text_length)
            call open_calculation_file (mcp_temporary_units(mcp_channel_index), output_file_prefix(1:prefix_length)//"."// &
                & number_text(1:2), "UNFORMATTED", &
                "UNKNOWN", open_status)
            if (open_status == 0) cycle
            do opened_channel_index = 30, mcp_channel_index - 1
                close(mcp_temporary_units(opened_channel_index))
            end do
            write (stderr_unit, *) "Error when opening the tmp files"
            stop
        end do

        do mcp_channel_index = 30, 32 + max_coulomb_multipole_rank
            write (mcp_temporary_units(mcp_channel_index)) "MCP", symmetry_block_index
        end do

        return
    end subroutine open_temporary_mcp_channels
end module open_temporary_mcp_channels_interface
