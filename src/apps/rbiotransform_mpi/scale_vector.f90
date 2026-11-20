!
!
!     ------------------------------------------------------------------
!     S C A L V E
!     ------------------------------------------------------------------
!
subroutine scale_vector(vector, factor, ndim)
    ! Module dependencies
    !
    ! CALCULATE SCALAR(FACTOR) TIMES VECTOR
    !************************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: ndim
    real(real64), intent(in) :: factor
    real(real64), intent(inout) :: vector(ndim)
    !-----------------------------------------------
    !
    vector(:ndim) = vector(:ndim)*factor
    !
    return
end subroutine scale_vector
