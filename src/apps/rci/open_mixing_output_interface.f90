module open_mixing_output_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine open_mixing_output (name, block_labels)
            implicit none (type, external)

            ! Arguments
            character (len = *), intent(in) :: name
            character (len = 8) :: block_labels(*)
        end subroutine open_mixing_output
    end interface
end module open_mixing_output_interface
