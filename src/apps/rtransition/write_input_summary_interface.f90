module write_input_summary_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_input_summary (name, inpci,ilbl)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: name(2)
            integer, intent(in) :: inpci
            integer :: ilbl
        end subroutine write_input_summary
    end interface
end module write_input_summary_interface
