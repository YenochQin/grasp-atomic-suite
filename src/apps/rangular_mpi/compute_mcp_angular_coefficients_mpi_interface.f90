module compute_mcp_angular_coefficients_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !cjb  myid, nprocs = NOT args
        subroutine compute_mcp_angular_coefficients_mpi (nb, is_restart, output_file_prefix)
            implicit none (type, external)

            ! Arguments
            integer  :: nb
            !     INTEGER  :: MYID
            !     INTEGER  :: NPROCS
            logical, intent(in) :: is_restart
            character(len=*), intent(in) :: output_file_prefix
        end subroutine compute_mcp_angular_coefficients_mpi
    end interface
end module compute_mcp_angular_coefficients_mpi_interface
