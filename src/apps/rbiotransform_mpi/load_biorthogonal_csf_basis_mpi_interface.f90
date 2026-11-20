module load_biorthogonal_csf_basis_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine load_biorthogonal_csf_basis_mpi (name, core_subshell_count,igg)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: name
            integer :: core_subshell_count,igg
        end subroutine load_biorthogonal_csf_basis_mpi
    end interface
end module load_biorthogonal_csf_basis_mpi_interface
