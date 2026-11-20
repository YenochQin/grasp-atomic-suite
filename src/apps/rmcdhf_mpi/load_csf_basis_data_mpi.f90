subroutine load_csf_basis_data_mpi(name, core_subshell_count, block_labels)
    ! Module dependencies
    use mpi_rmcdhf_csf_input, only: load_rmcdhf_descriptors_mpi
    use hamiltonian_symmetry_blocks, only: symmetry_block_count
    implicit none (type, external)

    ! Arguments
    character(len=*), intent(in) :: name
    integer, intent(out) :: core_subshell_count
    character(len=8), intent(in) :: block_labels(*)
    call load_rmcdhf_descriptors_mpi(name,core_subshell_count,block_labels(:symmetry_block_count))
end subroutine load_csf_basis_data_mpi
