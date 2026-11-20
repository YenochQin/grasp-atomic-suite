module try_allocate_sort_buffers_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine try_allocate_sort_buffers(n, irez)
            implicit none (type, external)

            ! Arguments
            integer, intent(in)  :: n
            integer, intent(out) :: irez
        end subroutine try_allocate_sort_buffers
    end interface
end module try_allocate_sort_buffers_interface
