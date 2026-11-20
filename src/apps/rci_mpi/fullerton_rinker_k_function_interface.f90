module fullerton_rinker_k_function_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function fullerton_rinker_k_function (x, n)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            real(real64), intent(in) :: x
            integer, intent(in) :: n
            !...This routine performs I/O.
        end function fullerton_rinker_k_function
    end interface
end module fullerton_rinker_k_function_interface
