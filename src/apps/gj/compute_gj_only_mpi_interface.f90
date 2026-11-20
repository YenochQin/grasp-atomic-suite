module compute_gj_only_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine compute_gj_only_mpi (gjc_diag, dgjc_diag)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use selected_level_indices, only: nvec
            implicit none (type, external)

            ! Arguments
            real(real64), intent(out) :: gjc_diag(nvec), dgjc_diag(nvec)
        end subroutine compute_gj_only_mpi
    end interface
end module compute_gj_only_mpi_interface
