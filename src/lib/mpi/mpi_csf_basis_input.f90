! Load the complete CSF basis on rank 0 and replicate it across ranks.
module mpi_csf_basis_input
    ! Module dependencies
    use mpi_csf_collectives, only: broadcast_descriptor_arrays, require_collective_success
    use mpi
    use mpi_state, only: mpi_rank
    use csf_full_basis_reader, only: load_csf_basis
    use grasp_parameters, only: max_orbitals
    use orbital_csf_basis, only: csf_count, orbital_count, csf_subshell_occupations, orbital_principal_n, orbital_kappa, &
        orbital_l, orbital_two_j, orbital_symmetry_labels
    use atomic_calculation_state, only: electron_count
    use csf_coupling_data, only: csf_subshell_quantum_codes, csf_angular_codes, csf_angular_code_count
    use transition_basis_blocks, only: transition_block_count, cumulative_csf_counts
    implicit none (type, external)

    ! Module interface
    private
    public :: load_csf_basis_mpi
contains

    subroutine load_csf_basis_mpi(descriptor_path,core_subshell_count)
        ! Arguments
        character(len=*), intent(in) :: descriptor_path
        integer, intent(out) :: core_subshell_count

        ! Local declarations
        integer :: descriptor_status, mpi_return_code, basis_counts(5)
        descriptor_status = 0
        if (mpi_rank == 0) call load_csf_basis(descriptor_path,core_subshell_count,descriptor_status)
        call require_collective_success(descriptor_status)
        if (mpi_rank == 0) basis_counts=[csf_count,orbital_count,electron_count,transition_block_count,core_subshell_count]
        call MPI_Bcast(basis_counts,5,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        csf_count=basis_counts(1)
        orbital_count=basis_counts(2)
        electron_count=basis_counts(3)
        transition_block_count=basis_counts(4)
        core_subshell_count=basis_counts(5)
        call MPI_Bcast(cumulative_csf_counts,size(cumulative_csf_counts),MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        ! Worker ranks allocate storage after receiving the basis dimensions.
        if (mpi_rank /= 0) then
            allocate(csf_subshell_occupations(max_orbitals,csf_count),stat=descriptor_status)
            if (descriptor_status == 0) allocate(csf_subshell_quantum_codes(max_orbitals,3,csf_count),stat=descriptor_status)
            if (descriptor_status == 0) allocate(csf_angular_codes(csf_angular_code_count,csf_count),stat=descriptor_status)
        end if
        call require_collective_success(descriptor_status)
        call MPI_Bcast(orbital_principal_n,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        call MPI_Bcast(orbital_kappa,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        call MPI_Bcast(orbital_l,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        call MPI_Bcast(orbital_two_j,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        call MPI_Bcast(orbital_symmetry_labels,2*orbital_count,MPI_CHARACTER,0,MPI_COMM_WORLD,mpi_return_code)
        call broadcast_descriptor_arrays(csf_subshell_occupations,csf_subshell_quantum_codes,csf_angular_codes)
    end subroutine load_csf_basis_mpi
end module mpi_csf_basis_input
