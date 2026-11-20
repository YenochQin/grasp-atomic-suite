module write_input_summary_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_input_summary
            implicit none (type, external)
            !...This routine performs I/O.
        end subroutine write_input_summary
    end interface
end module write_input_summary_interface
