module load_merged_csf_basis_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine load_merged_csf_basis (name, core_subshell_count,igg)
            implicit none (type, external)

            ! Arguments
            character (len = 128), intent(in) :: name
            integer :: core_subshell_count,igg
        end subroutine load_merged_csf_basis
    end interface
end module load_merged_csf_basis_interface
