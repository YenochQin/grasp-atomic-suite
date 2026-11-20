module open_summary_file_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine open_summary_file(default_file_path)
            implicit none (type, external)

            ! Arguments
            character(len=*), intent(in) :: default_file_path
        end subroutine open_summary_file
    end interface
end module open_summary_file_interface
