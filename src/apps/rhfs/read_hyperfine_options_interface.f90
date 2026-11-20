module read_hyperfine_options_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_hyperfine_options (name)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: name
        end subroutine read_hyperfine_options
    end interface
end module read_hyperfine_options_interface
