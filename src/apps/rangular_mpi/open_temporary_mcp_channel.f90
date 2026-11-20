!***********************************************************************
!                                                                      *
subroutine open_temporary_mcp_channel(nb, kmax, filehead)
    ! Module dependencies
    use mcp_generation_options, only: mcp_temporary_units
    !***********************************************************************
    use integer_format, only: convrt
    use file_open, only: openfl
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    implicit none (type, external)

    ! Arguments
    integer , intent(in) :: nb
    integer , intent(in) :: kmax
    character , intent(in) :: filehead*(*)

    ! Local declarations
    integer :: k, text_length, open_status, i, lng
    character :: number_text*2
    !-----------------------------------------------
    !
    ! All files  filehead.XX  are UNFORMATTED;
    !
    lng = len_trim(filehead)
    do k = 30, 32 + kmax
        call convrt (k, number_text, text_length)
        call openfl (mcp_temporary_units(k), filehead(1:lng)//"."//number_text(1:2), "UNFORMATTED", &
            "UNKNOWN", open_status)
        if (open_status == 0) cycle
        do i = 30, k - 1
            close(mcp_temporary_units(i))
        end do
        write (stderr_unit, *) "Error when opening the tmp files"
        stop
    end do

    do k = 30, 32 + kmax
        write (mcp_temporary_units(k)) "MCP", nb
    end do

    return
end subroutine open_temporary_mcp_channel
