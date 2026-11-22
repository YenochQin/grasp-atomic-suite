! Load RMCDHF CSF descriptors on rank 0 and distribute orbital and occupation data.
module mpi_rmcdhf_csf_input
    ! Module dependencies
    use mpi_csf_collectives, only: broadcast_descriptor_arrays, require_collective_success
    use mpi
    use mpi_state, only: mpi_rank
    use rmcdhf_csf_input, only: load_rmcdhf_descriptors
    use grasp_parameters, only: max_orbitals
    use orbital_csf_basis, only: csf_count, orbital_count, csf_subshell_occupations, orbital_principal_n, orbital_kappa, &
        orbital_l, orbital_two_j, orbital_symmetry_labels
    use atomic_calculation_state, only: electron_count
    use closed_core_subshells, only: closed_core_subshell_count
    use hamiltonian_symmetry_blocks, only: symmetry_block_count, csf_counts_by_block
    use level_symmetry_metadata, only: packed_block_symmetries, max_symmetry_blocks
    implicit none (type, external)

    ! Module interface
    private
    public :: load_rmcdhf_descriptors_mpi
contains

    subroutine load_rmcdhf_descriptors_mpi(descriptor_path,core_subshell_count,block_labels,expected_basis_id)
        ! Arguments
        character(len=*), intent(in) :: descriptor_path
        integer, intent(out) :: core_subshell_count
        character(len=8), intent(in) :: block_labels(:)
        character(len=64), intent(in), optional :: expected_basis_id

        ! Local declarations
        integer :: descriptor_status, mpi_return_code, basis_counts(5), root_csf_counts_by_block(max_symmetry_blocks)
        character(len=8) :: root_block_labels(max_symmetry_blocks)
        descriptor_status = 0
        if (mpi_rank == 0) call load_rmcdhf_descriptors(descriptor_path,core_subshell_count,block_labels,descriptor_status, &
            expected_basis_id)
        call require_collective_success(descriptor_status)
        if (mpi_rank == 0) then
            basis_counts = [csf_count,orbital_count,electron_count,symmetry_block_count,core_subshell_count]
            root_csf_counts_by_block = 0
            root_csf_counts_by_block(:symmetry_block_count) = csf_counts_by_block(1:symmetry_block_count)
            root_block_labels = ""
            root_block_labels(:symmetry_block_count) = block_labels
        end if
        ! Check the descriptor counts and block labels against the state already initialized on each rank.
        call MPI_Bcast(basis_counts,5,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        if (any(basis_counts /= [csf_count,orbital_count,electron_count,symmetry_block_count,closed_core_subshell_count])) &
            descriptor_status=1
        call require_collective_success(descriptor_status)
        call MPI_Bcast(root_csf_counts_by_block,max_symmetry_blocks,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        call MPI_Bcast(root_block_labels,8*max_symmetry_blocks,MPI_CHARACTER,0,MPI_COMM_WORLD,mpi_return_code)
        if (any(root_csf_counts_by_block(:symmetry_block_count) /= csf_counts_by_block(1:symmetry_block_count))) descriptor_status=1
        if (any(root_block_labels(:symmetry_block_count) /= block_labels)) descriptor_status=1
        call require_collective_success(descriptor_status)
        core_subshell_count=basis_counts(5)
        if (mpi_rank /= 0) allocate(csf_subshell_occupations(max_orbitals,csf_count),stat=descriptor_status)
        call require_collective_success(descriptor_status)
        call MPI_Bcast(orbital_principal_n,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        call MPI_Bcast(orbital_kappa,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        call MPI_Bcast(orbital_l,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        call MPI_Bcast(orbital_two_j,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        call MPI_Bcast(orbital_symmetry_labels,2*orbital_count,MPI_CHARACTER,0,MPI_COMM_WORLD,mpi_return_code)
        call MPI_Bcast(packed_block_symmetries,max_symmetry_blocks,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        call broadcast_descriptor_arrays(csf_subshell_occupations)
    end subroutine load_rmcdhf_descriptors_mpi
end module mpi_rmcdhf_csf_input
