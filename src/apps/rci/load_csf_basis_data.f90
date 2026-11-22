! Load descriptor data for the application's CSF basis and check its block metadata.
module load_csf_basis_data_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: load_csf_basis_data
contains

    subroutine load_csf_basis_data(name,core_subshell_count,input_symmetry_block_count,block_labels)
        ! Module dependencies
        use rci_csf_block_reader, only: open_csf_blocks
        implicit none (type, external)

        ! Arguments
        character(len=*) :: name
        integer :: core_subshell_count
        integer, intent(in) :: input_symmetry_block_count
        character(len=8) :: block_labels(*)

        ! Local declarations
        integer :: status
        call open_csf_blocks(name,core_subshell_count,block_labels(:input_symmetry_block_count),status)
        if (status /= 0) error stop "RCI descriptor header failed"
    end subroutine load_csf_basis_data
end module load_csf_basis_data_interface
