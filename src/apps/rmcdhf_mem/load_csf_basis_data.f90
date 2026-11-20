module load_csf_basis_data_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: load_csf_basis_data
contains

    subroutine load_csf_basis_data(name, core_subshell_count, block_labels)
        ! Module dependencies
        use rmcdhf_csf_input, only: load_rmcdhf_descriptors
        use hamiltonian_symmetry_blocks, only: symmetry_block_count
        implicit none (type, external)

        ! Arguments
        character(len=*), intent(in) :: name
        integer :: core_subshell_count
        character(len=8) :: block_labels(*)

        ! Local declarations
        integer :: status
        call load_rmcdhf_descriptors(name,core_subshell_count,block_labels(:symmetry_block_count),status)
        if (status /= 0) error stop "RMCDHF descriptor input failed"
    end subroutine load_csf_basis_data
end module load_csf_basis_data_interface
