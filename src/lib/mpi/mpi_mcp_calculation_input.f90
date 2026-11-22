! Install the indexed MCP metadata shared by the MPI MCDHF solvers.
module mpi_mcp_calculation_input
    use mpi_indexed_mcp, only: indexed_mcp, open_indexed_mcp_reader
    use array_memory_management, only: alloc
    use atomic_calculation_state, only: electron_count
    use first_order_csf_cutoffs, only: first_order_csf_cutoff
    use hamiltonian_symmetry_blocks, only: csf_counts_by_block, symmetry_block_count
    use mcp_file_metadata, only: max_mcp_file_multipole_rank
    use orbital_csf_basis, only: csf_count, orbital_count
    use rmcdhf_mcp_options, only: mcp_diagonal_only, mcp_first_order_perturbation
    implicit none (type, external)
    private
    public :: indexed_input, open_calculation_mcp

    type(indexed_mcp) :: indexed_input
contains
    subroutine open_calculation_mcp(core_subshell_count, maximum_blocks, block_labels, prefix)
        integer, intent(out) :: core_subshell_count
        integer, intent(in) :: maximum_blocks
        character(len=8), intent(out) :: block_labels(maximum_blocks)
        character(*), intent(in) :: prefix

        call open_indexed_mcp_reader(indexed_input,prefix,maximum_blocks)
        core_subshell_count = indexed_input%core_count
        symmetry_block_count = indexed_input%block_count
        max_mcp_file_multipole_rank = indexed_input%max_rank
        electron_count = indexed_input%electron_count
        csf_count = indexed_input%csf_count
        orbital_count = indexed_input%orbital_count
        mcp_diagonal_only = btest(indexed_input%interaction_flags,0)
        mcp_first_order_perturbation = btest(indexed_input%interaction_flags,1)
        first_order_csf_cutoff = int(indexed_input%cutoffs(1))
        call alloc(csf_counts_by_block,0,symmetry_block_count,"csf_counts_by_block","open_calculation_mcp")
        csf_counts_by_block(0) = 0
        csf_counts_by_block(1:symmetry_block_count) = int(indexed_input%block_sizes)
        block_labels = ""
        block_labels(1:symmetry_block_count) = indexed_input%block_labels
    end subroutine open_calculation_mcp
end module mpi_mcp_calculation_input
