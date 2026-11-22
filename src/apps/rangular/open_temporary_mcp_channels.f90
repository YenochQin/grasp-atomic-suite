module open_temporary_mcp_channels_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_temporary_mcp_channels
contains

    ! Modified by Gediminas Gaigalas: Feb 2017
    !                             1) for new spin-angular integration,
    !                             2) for sorting in the memory.
    subroutine open_temporary_mcp_channels(nb, k, filehead)
        ! Module dependencies
        use mcp_generation_options, only: mcp_temporary_units
        use integer_format, only: format_integer
        use file_open, only: open_calculation_file
        use standard_io_units, only: stderr_unit
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: nb
        integer, intent(in) :: k
        character, intent(in) :: filehead*(*)

        ! Local declarations
        integer   :: text_length, open_status, lng
        character :: number_text*2
        ! All files  filehead.XX  are UNFORMATTED;
        lng = len_trim(filehead)
        call format_integer (k, number_text, text_length)
        call open_calculation_file (mcp_temporary_units(k), filehead(1:lng)//"."//number_text(1:2), "UNFORMATTED", &
            "UNKNOWN", open_status)
        if (open_status /= 0) then
            write (stderr_unit, *) "Error when opening the tmp files"
            stop
        end if
        write (mcp_temporary_units(k)) "MCP", nb
        return
    end subroutine open_temporary_mcp_channels
end module open_temporary_mcp_channels_interface
