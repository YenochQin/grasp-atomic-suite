module configure_debug_output_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine configure_debug_output(default_file_path)
            implicit none (type, external)

            ! Arguments
            character(len=*), intent(in) :: default_file_path
        end subroutine configure_debug_output
    end interface
end module configure_debug_output_interface
