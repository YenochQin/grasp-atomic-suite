module electron_screening_factor_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function electron_screening_factor (j)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: j
        end function electron_screening_factor
    end interface
end module electron_screening_factor_interface
