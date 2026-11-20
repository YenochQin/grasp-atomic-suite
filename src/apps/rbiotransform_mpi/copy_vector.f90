!
!     ------------------------------------------------------------------
!     C O P V E C
!     ------------------------------------------------------------------
!
subroutine copy_vector(from, to, ndim)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: ndim
    real(real64), intent(in)  :: from(ndim)
    real(real64), intent(out) :: to(ndim)
    !
    to(:ndim) = from(:ndim)
    !
    return
end subroutine copy_vector
