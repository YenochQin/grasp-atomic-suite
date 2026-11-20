module configure_debug_output_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine configure_debug_output_mpi (debug, fullname)
            implicit none (type, external)

            ! Arguments
            logical :: debug
            character(len=*) :: fullname
        end subroutine configure_debug_output_mpi
    end interface
end module configure_debug_output_mpi_interface
