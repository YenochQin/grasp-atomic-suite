module vector_inner_product_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function vector_inner_product (a, b, ndim)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            real(real64), intent(in) :: a(*)
            real(real64), intent(in) :: b(*)
            integer, intent(in) :: ndim
        end function vector_inner_product
    end interface
end module vector_inner_product_interface
