module load_csf_basis_data_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine load_csf_basis_data (name, core_subshell_count, block_labels)
            implicit none (type, external)

            ! Arguments
            character (len = *), intent(in) :: name
            integer :: core_subshell_count
            character (len = 8) :: block_labels(*)
        end subroutine load_csf_basis_data
    end interface
end module load_csf_basis_data_interface
