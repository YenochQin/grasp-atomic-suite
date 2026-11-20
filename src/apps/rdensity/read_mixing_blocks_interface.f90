module read_mixing_blocks_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_mixing_blocks (state_name, ci_input_mode)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: state_name
            integer, intent(in) :: ci_input_mode
        end subroutine read_mixing_blocks
    end interface
end module read_mixing_blocks_interface
