module configure_debug_output_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine configure_debug_output_mpi (debug_file_path)
            implicit none (type, external)

            ! Arguments
            character (len = *), intent(in) :: debug_file_path
        end subroutine configure_debug_output_mpi
    end interface
end module configure_debug_output_mpi_interface
