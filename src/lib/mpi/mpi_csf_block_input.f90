! Keep the RCI descriptor stream on rank 0 and distribute one symmetry block at a time.
module mpi_csf_block_input
    ! Module dependencies
    use rci_csf_block_reader, only: open_csf_blocks, load_csf_block
    use mpi_csf_collectives, only: require_collective_success, broadcast_descriptor_arrays
    use mpi_state, only: mpi_rank, mpi_bcast, mpi_integer, mpi_character, mpi_comm_world
    use orbital_csf_basis, only: csf_count, orbital_count, orbital_principal_n, orbital_kappa, orbital_l, orbital_two_j, &
        orbital_symmetry_labels, csf_subshell_occupations
    use csf_coupling_data, only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
    use atomic_calculation_state, only: electron_count
    use hamiltonian_symmetry_blocks, only: symmetry_block_count, csf_counts_by_block
    implicit none (type, external)

    ! Module interface
    private
    public :: open_csf_blocks_mpi, load_csf_block_mpi
contains
    ! Distribute orbital and block metadata before any CSF rows are loaded.
    subroutine open_csf_blocks_mpi(descriptor_path,core_subshell_count,block_labels)
        ! Arguments
        character(len=*), intent(in) :: descriptor_path
        integer, intent(out) :: core_subshell_count
        character(len=8), intent(out) :: block_labels(:)

        ! Local declarations
        integer :: descriptor_status, mpi_return_code, basis_counts(5)
        descriptor_status=0
        if (mpi_rank == 0) call open_csf_blocks(descriptor_path,core_subshell_count,block_labels,descriptor_status)
        call require_collective_success(descriptor_status)
        if (mpi_rank == 0) basis_counts=[csf_count,orbital_count,core_subshell_count,electron_count,symmetry_block_count]
        call mpi_bcast(basis_counts,5,mpi_integer,0,mpi_comm_world,mpi_return_code)
        call require_collective_success(mpi_return_code)
        csf_count=basis_counts(1)
        orbital_count=basis_counts(2)
        core_subshell_count=basis_counts(3)
        electron_count=basis_counts(4)
        symmetry_block_count=basis_counts(5)
        call require_collective_success(merge(0,1,size(block_labels) >= symmetry_block_count))
        if (mpi_rank /= 0) then
            if (associated(csf_counts_by_block)) deallocate(csf_counts_by_block)
            allocate(csf_counts_by_block(0:symmetry_block_count),stat=descriptor_status)
            block_labels=""
        end if
        call require_collective_success(descriptor_status)
        call mpi_bcast(csf_counts_by_block,symmetry_block_count+1,mpi_integer,0,mpi_comm_world,mpi_return_code)
        call mpi_bcast(block_labels,8*symmetry_block_count,mpi_character,0,mpi_comm_world,mpi_return_code)
        call mpi_bcast(orbital_principal_n,orbital_count,mpi_integer,0,mpi_comm_world,mpi_return_code)
        call mpi_bcast(orbital_kappa,orbital_count,mpi_integer,0,mpi_comm_world,mpi_return_code)
        call mpi_bcast(orbital_l,orbital_count,mpi_integer,0,mpi_comm_world,mpi_return_code)
        call mpi_bcast(orbital_two_j,orbital_count,mpi_integer,0,mpi_comm_world,mpi_return_code)
        call mpi_bcast(orbital_symmetry_labels,2*orbital_count,mpi_character,0,mpi_comm_world,mpi_return_code)
        call require_collective_success(mpi_return_code)
    end subroutine open_csf_blocks_mpi

    subroutine load_csf_block_mpi(block_index)
        ! Arguments
        integer, intent(in) :: block_index

        ! Local declarations
        integer :: descriptor_status
        descriptor_status=0
        if (mpi_rank == 0) call load_csf_block(block_index,descriptor_status)
        call require_collective_success(descriptor_status)
        call broadcast_descriptor_arrays(csf_subshell_occupations,csf_subshell_quantum_codes,csf_coupling_quantum_codes)
    end subroutine load_csf_block_mpi
end module mpi_csf_block_input
