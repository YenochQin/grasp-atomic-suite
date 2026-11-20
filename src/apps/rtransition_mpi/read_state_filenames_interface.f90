module read_state_filenames_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_state_filenames (name)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(out) :: name(2)
            !...This routine performs I/O.
        end subroutine read_state_filenames
    end interface
end module read_state_filenames_interface
