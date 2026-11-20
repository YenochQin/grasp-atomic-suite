module read_eol_average_options_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_eol_average_options_mpi (idblk)
            implicit none (type, external)

            ! Arguments
            character (len = 8) :: idblk(*)
        end subroutine read_eol_average_options_mpi
    end interface
end module read_eol_average_options_mpi_interface
