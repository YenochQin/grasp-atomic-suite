! Open this rank's MCP output channels and record the MPI layout in their headers.
module open_mcp_outputs_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_mcp_outputs_mpi
contains

    !  Set interaction options and write the MPI MCP file headers.
    !
    ! Written by Xinghong He. Last revision: 30 Jun 1998
    !
    ! cjb mpi_C
    subroutine open_mcp_outputs_mpi(core_subshell_count, block_labels, filehead, descriptors, indexed_output)
        ! Module dependencies
        use mpi_indexed_mcp, only: indexed_mcp, open_indexed_mcp_writer
        use hamiltonian_symmetry_blocks, only: symmetry_block_count, csf_counts_by_block
        use csf_descriptor_reader, only: descriptor_source
        use mpi_csf_collectives, only: require_collective_success
        use csf_block_cutoffs, only: active_csf_cutoffs
        use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
        use orbital_csf_basis, only: csf_count, orbital_count
        use atomic_calculation_state, only: electron_count
        ! cjb mpi_C
        use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_character, mpi_integer, mpi_logical, mpi_rank

        ! Procedure interfaces
        use connect_mcp_files_interface, only: configure_mcp_channels
        use read_calculation_options_interface, only: read_calculation_options
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: core_subshell_count
        character , intent(in) :: filehead*(*)
        character , intent(in) :: block_labels(*)*8
        type(descriptor_source), intent(in) :: descriptors
        type(indexed_mcp), intent(out) :: indexed_output

        !  Local declarations
        ! cjb mpi_C
        integer :: input_status
        character(len=64) :: basis_id

        input_status = 0
        if (mpi_rank == 0) call read_calculation_options(descriptors, input_status)
        call require_collective_success(input_status)

        call mpi_bcast (active_csf_cutoffs, size(active_csf_cutoffs), mpi_integer, 0, mpi_comm_world, mpi_error_code)
        call mpi_bcast (diagonal_only,  1,mpi_logical,0,mpi_comm_world,mpi_error_code)
        call mpi_bcast (use_first_order_restriction,1,mpi_logical,0,mpi_comm_world,mpi_error_code)

        call configure_mcp_channels()
        basis_id = descriptors%basis_id
        call mpi_bcast(basis_id,len(basis_id),mpi_character,0,mpi_comm_world,mpi_error_code)
        call open_indexed_mcp_writer(indexed_output,filehead,core_subshell_count,max_coulomb_multipole_rank, &
            electron_count,orbital_count,csf_counts_by_block(1:symmetry_block_count),block_labels(:symmetry_block_count), &
            active_csf_cutoffs(:symmetry_block_count),merge(1,0,diagonal_only)+merge(2,0,use_first_order_restriction),basis_id)
    end subroutine open_mcp_outputs_mpi
end module open_mcp_outputs_mpi_interface
