! Manage sequential CSF descriptor blocks for RCI matrix construction.
module rci_csf_block_reader
    ! Module dependencies
    use csf_descriptor_reader, only: descriptor_source, open_descriptor_source, close_descriptor_source, &
        install_descriptor_header, read_descriptor_block, read_descriptor_zero_order_counts
    use rmcdhf_csf_input, only: descriptor_path
    use hamiltonian_symmetry_blocks, only: csf_counts_by_block
    use orbital_csf_basis, only: csf_subshell_occupations
    use csf_coupling_data, only: csf_subshell_quantum_codes, csf_angular_codes
    implicit none (type, external)

    ! Module interface
    private

    ! Module state
    ! The stream persists between calls while RCI constructs one symmetry block at a time.
    type(descriptor_source), save :: block_reader

    ! Module interface
    public :: open_csf_blocks, load_csf_block
    public :: read_csf_zero_order_counts
contains
    subroutine open_csf_blocks(state_name,core_subshell_count,block_labels,status)
        ! Arguments
        character(len=*), intent(in) :: state_name
        integer, intent(out) :: core_subshell_count,status
        character(len=8), intent(out) :: block_labels(:)
        call open_descriptor_source(block_reader,descriptor_path(state_name),size(block_labels),status=status,spin_angular=.true.)
        if (status /= 0) return
        call install_descriptor_header(block_reader,core_subshell_count,block_labels,status)
        if (status == 0) then
            deallocate(csf_counts_by_block)
            allocate(csf_counts_by_block(0:block_reader%block_count),stat=status)
            if (status == 0) then
                csf_counts_by_block(0)=0
                csf_counts_by_block(1:)=block_reader%block_lengths
            end if
        end if
        if (status /= 0) call close_descriptor_source(block_reader)
    end subroutine open_csf_blocks

    subroutine read_csf_zero_order_counts(zero_order_counts,status)
        ! Arguments
        integer, intent(out) :: zero_order_counts(:), status
        call read_descriptor_zero_order_counts(block_reader,zero_order_counts,status)
    end subroutine read_csf_zero_order_counts

    ! Replace the active block arrays and close the reader at the last block or on failure.
    subroutine load_csf_block(block_index,status)
        ! Arguments
        integer, intent(in) :: block_index
        integer, intent(out) :: status
        call read_descriptor_block(block_reader,block_index,csf_subshell_occupations, &
            csf_subshell_quantum_codes,csf_angular_codes,status)
        if (status /= 0 .or. block_index == block_reader%block_count) call close_descriptor_source(block_reader)
    end subroutine load_csf_block
end module rci_csf_block_reader
