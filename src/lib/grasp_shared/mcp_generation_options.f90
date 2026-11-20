module mcp_generation_options
    implicit none (type, external)

    ! Module state
    integer :: max_coulomb_multipole_rank
    logical :: diagonal_only, use_first_order_restriction
    ! File suffixes index the output and per-block temporary connections.
    integer, allocatable :: mcp_output_units(:), mcp_temporary_units(:)
end module mcp_generation_options
