module query_total_system_memory_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function query_total_system_memory ()
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)
            !...This routine performs I/O.
        end function query_total_system_memory
    end interface
end module query_total_system_memory_interface
