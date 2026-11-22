! Read the rank-local indexed MCP metadata before loading the basis.
module connect_mcp_files_interface
    use mpi_mcp_calculation_input, only: connect_mcp_files => open_calculation_mcp
    implicit none (type, external)
    private
    public :: connect_mcp_files
end module connect_mcp_files_interface
