!
!     ------------------------------------------------------------------
!     I N P R O D
!     ------------------------------------------------------------------
!
real(real64) function vector_inner_product (a, b, ndim)
    ! Module dependencies
    !      CALCULATE SCALAR PRODUCT BETWEEN TO VECTORS A,B
    !************************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: ndim
    real(real64), intent(in) :: a(*)
    real(real64), intent(in) :: b(*)

    ! Local declarations
    integer :: i
    !-----------------------------------------------
    !
    vector_inner_product = dot_product(a(:ndim),b(:ndim))
    !
    return
end function vector_inner_product
