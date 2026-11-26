! Compile the unmodified production MATRIX body; only assembly/eigensolution
! and memory-record input are controlled seams. BLAS is the linked real library.
module ci_fixture
    use eigenstate_expansions, only: csf_mixing_coefficients, shifted_eigenvalues
    implicit none (type, external)
    private
    public :: shifted_eigenvalues, csf_mixing_coefficients, solves, eigenpairs
    integer :: solves = 0
contains
    subroutine eigenpairs()
        use, intrinsic :: iso_fortran_env, only: real64
        solves = solves + 1
        shifted_eigenvalues = [-2.0_real64, -1.0_real64]
        if (solves == 1) then
            csf_mixing_coefficients(3:6) = [1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64]
        else
            csf_mixing_coefficients(3:6) = [0.6_real64, 0.8_real64, -0.8_real64, 0.6_real64]
        end if
    end subroutine eigenpairs
end module ci_fixture

! Supply module-level substitutes without linking the corresponding app implementations.
#ifdef CI_MPI
module diagonalize_hamiltonian_block_mpi_interface
#else
module diagonalize_hamiltonian_block_interface
#endif
    implicit none (type, external)
    private
#ifdef CI_MPI
    public :: diagonalize_hamiltonian_block_mpi
#else
    public :: diagonalize_hamiltonian_block
#endif
contains
#ifdef CI_MPI
    subroutine diagonalize_hamiltonian_block_mpi(dvdfirst, print_details, symmetry_block_index, &
        ncfpat, ncminpat, nevecpat, total_csf_count)
#else
    subroutine diagonalize_hamiltonian_block(dvdfirst, print_details, symmetry_block_index, &
        ncfpat, ncminpat, nevecpat, total_csf_count)
#endif
        use ci_fixture, only: eigenpairs
        implicit none (type, external)
        logical, intent(in) :: dvdfirst
        logical :: print_details
        integer, intent(in) :: symmetry_block_index, ncfpat, ncminpat, nevecpat
        integer :: total_csf_count
        if (symmetry_block_index /= 1 .or. nevecpat /= 2 .or. total_csf_count /= 2) error stop "eigensolver seam arguments"
        call eigenpairs()
#ifdef CI_MPI
    end subroutine diagonalize_hamiltonian_block_mpi
end module diagonalize_hamiltonian_block_mpi_interface
#else
    end subroutine diagonalize_hamiltonian_block
end module diagonalize_hamiltonian_block_interface
#endif

module assemble_hamiltonian_block_interface
    implicit none (type, external)
    private
    public :: assemble_hamiltonian_block
contains
#ifdef CI_MEM
    subroutine assemble_hamiltonian_block(symmetry_block_index, mpi_rank, mpi_process_count, &
        total_contribution_count, total_secondary_contribution_count)
#else
    subroutine assemble_hamiltonian_block(symmetry_block_index, mpi_rank, mpi_process_count)
#endif
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use sparse_hamiltonian_storage, only: hamiltonian_values, matrix_partner_indices, stored_segment_end_offsets
        implicit none (type, external)
        integer, intent(in) :: symmetry_block_index, mpi_rank, mpi_process_count
#ifdef CI_MEM
        integer(int64) :: total_contribution_count(20), total_secondary_contribution_count(20)
#endif
        integer :: i
        if (symmetry_block_index /= 1) error stop "Hamiltonian seam block"
        if (.not. allocated(hamiltonian_values) .or. .not. allocated(matrix_partner_indices)) &
            error stop "assembly needs allocated arrays even on empty ranks"
        if (size(hamiltonian_values) /= count([(mod(i - 1, mpi_process_count) == mpi_rank, i=1,2)]) .or. &
            size(matrix_partner_indices) /= size(hamiltonian_values)) error stop "wrong rank-local storage size"
        if (lbound(stored_segment_end_offsets,1) /= 0 .or. stored_segment_end_offsets(0) /= 0_int64) &
            error stop "missing zero-based offset sentinel"
        hamiltonian_values = 0.0_real64
    end subroutine assemble_hamiltonian_block
end module assemble_hamiltonian_block_interface

#ifdef CI_MEM
module read_mcp_coefficients_from_memory_interface
    implicit none (type, external)
    private
    public :: read_mcp_coefficients_from_memory
contains
    subroutine read_mcp_coefficients_from_memory(jblock)
        use sparse_hamiltonian_storage, only: matrix_partner_indices, stored_segment_end_offsets
#ifdef CI_MPI
        use mpi_state, only: mpi_rank, mpi_process_count
#else
        use serial_mpi_context, only: mpi_rank, mpi_process_count
#endif
        implicit none (type, external)
        integer, intent(in) :: jblock
        integer :: i, k
        if (jblock /= 1) error stop "memory seam block"
        k = 0
        do i = mpi_rank + 1, 2, mpi_process_count
            k = k + 1
            stored_segment_end_offsets(i) = k
            matrix_partner_indices(k) = i
        end do
    end subroutine read_mcp_coefficients_from_memory
end module read_mcp_coefficients_from_memory_interface

module load_mcp_records_into_memory_interface
    implicit none (type, external)
    private
    public :: load_mcp_records_into_memory
contains
    subroutine load_mcp_records_into_memory(nfile, total_secondary_contribution_count, lab, ncontr)
        use, intrinsic :: iso_fortran_env, only: int64
        implicit none (type, external)
        integer, intent(in) :: nfile
        integer(int64), intent(in) :: total_secondary_contribution_count
        integer, intent(out) :: lab, ncontr
        if (nfile < 31 .or. nfile > 33 .or. total_secondary_contribution_count /= 1) &
            error stop "unexpected empty block seam"
        lab = 0
        ncontr = 0
    end subroutine load_mcp_records_into_memory
end module load_mcp_records_into_memory_interface
#endif
