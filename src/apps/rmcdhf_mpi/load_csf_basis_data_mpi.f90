! Load the MCDHF descriptor basis on rank 0 and distribute it after checking the MCP metadata.
module load_csf_basis_data_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: load_csf_basis_data_mpi
contains

    subroutine load_csf_basis_data_mpi(name, core_subshell_count, block_labels)
        ! Module dependencies
        use mpi_rmcdhf_csf_input, only: load_rmcdhf_descriptors_mpi
        use hamiltonian_symmetry_blocks, only: symmetry_block_count
        implicit none (type, external)

        ! Arguments
        character(len=*), intent(in) :: name
        integer :: core_subshell_count
        character(len=8) :: block_labels(*)
        call load_rmcdhf_descriptors_mpi(name,core_subshell_count,block_labels(:symmetry_block_count))
    end subroutine load_csf_basis_data_mpi
end module load_csf_basis_data_mpi_interface
