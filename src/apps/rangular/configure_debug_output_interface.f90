module configure_debug_output_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine configure_debug_output (debug, file_path)
            implicit none (type, external)

            ! Arguments
            logical :: debug
            character(len=*) :: file_path
        end subroutine configure_debug_output
    end interface
end module configure_debug_output_interface
