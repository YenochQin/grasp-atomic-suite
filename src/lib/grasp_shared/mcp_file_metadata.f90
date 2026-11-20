module mcp_file_metadata
    implicit none (type, external)

    ! Module state
    integer :: max_mcp_file_multipole_rank
    ! Indexed by the numeric MCP filename suffix, independent of the runtime unit.
    integer, allocatable :: mcp_input_units(:)
end module mcp_file_metadata
