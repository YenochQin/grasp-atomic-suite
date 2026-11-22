! Configure the rank-local generation channels from the orbital basis.
module connect_mcp_files_interface
    use mcp_generation_options, only: mcp_temporary_units, max_coulomb_multipole_rank
    use orbital_csf_basis, only: orbital_count, orbital_two_j
    use mpi_io_errors, only: fail_input_io_mpi
    implicit none (type, external)
    private
    public :: configure_mcp_channels
contains
    subroutine configure_mcp_channels()
        max_coulomb_multipole_rank = maxval(orbital_two_j(:orbital_count))
        if (max_coulomb_multipole_rank > 67) &
            call fail_input_io_mpi(0,-1,"MCP temporary channel exceeds two-digit suffix")
        allocate(mcp_temporary_units(30:32+max_coulomb_multipole_rank))
    end subroutine configure_mcp_channels
end module connect_mcp_files_interface
