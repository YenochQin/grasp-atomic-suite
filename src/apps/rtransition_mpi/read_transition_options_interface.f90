module read_transition_options_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_transition_options (name)
            implicit none (type, external)

            ! Arguments
            character (len = 128) :: name(2)
        end subroutine read_transition_options
    end interface
end module read_transition_options_interface
