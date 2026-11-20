module fullerton_rinker_l_function_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function fullerton_rinker_l_function (x, k)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            real(real64), intent(in) :: x
            integer, intent(in) :: k
            !...This routine performs I/O.
        end function fullerton_rinker_l_function
    end interface
end module fullerton_rinker_l_function_interface
