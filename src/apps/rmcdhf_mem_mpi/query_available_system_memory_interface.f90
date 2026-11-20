module query_available_system_memory_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !   Written by  G. Gaigalas               Vilnius, September 2021
        real(real64) function query_available_system_memory ()
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)
        end function query_available_system_memory
    end interface
end module query_available_system_memory_interface
