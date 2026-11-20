module load_csf_basis_data_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine load_csf_basis_data_mpi (name, core_subshell_count, block_labels)
            implicit none (type, external)

            ! Arguments
            character (len = *), intent(in) :: name
            integer :: core_subshell_count
            character (len = 8) :: block_labels(*)
            !VAST...Dummy argument IDBLK is not referenced in this routine.
        end subroutine load_csf_basis_data_mpi
    end interface
end module load_csf_basis_data_mpi_interface
