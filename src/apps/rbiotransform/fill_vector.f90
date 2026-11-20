!
!     ------------------------------------------------------------------
!     S E T V E C
!     ------------------------------------------------------------------
!
subroutine fill_vector(vector, value, ndim)
    ! Module dependencies
    !
    ! VECTOR (*) = VALUE
    !
    !************************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: ndim
    real(real64), intent(in) :: value
    real(real64), intent(out) :: vector(ndim)

    ! Local declarations
    integer :: i
    !-----------------------------------------------
    !
    vector(:ndim) = value
    !
    return
end subroutine fill_vector
