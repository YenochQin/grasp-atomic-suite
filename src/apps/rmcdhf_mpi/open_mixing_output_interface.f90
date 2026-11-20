module open_mixing_output_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine open_mixing_output (name)
            implicit none (type, external)

            ! Arguments
            character (len = *), intent(in) :: name
        end subroutine open_mixing_output
    end interface
end module open_mixing_output_interface
