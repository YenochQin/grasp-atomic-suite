module hyperfine_one_electron_element_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine hyperfine_one_electron_element (i1, tensor_rank, i2, apart, gjpart, dgjpart)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: i1
            integer, intent(in) :: tensor_rank
            integer, intent(in) :: i2
            real(real64), intent(out) :: apart
            real(real64), intent(out) :: gjpart
            real(real64), intent(out) :: dgjpart
        end subroutine hyperfine_one_electron_element
    end interface
end module hyperfine_one_electron_element_interface
