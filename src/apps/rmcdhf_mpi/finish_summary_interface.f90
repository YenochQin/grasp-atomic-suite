module finish_summary_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine finish_summary
            implicit none (type, external)
            !...This routine performs I/O.
        end subroutine finish_summary
    end interface
end module finish_summary_interface
