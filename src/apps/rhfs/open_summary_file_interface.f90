module open_summary_file_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine open_summary_file (name, nci)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: name
            integer, intent(in) :: nci
        end subroutine open_summary_file
    end interface
end module open_summary_file_interface
