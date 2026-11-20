!
!     ------------------------------------------------------------------
!     T R P M A T
!     ------------------------------------------------------------------
!
subroutine transpose_matrix(xin, nrow, ncol, xout)
    ! Module dependencies
    !
    ! XOUT(I,J) = XIN(J,I)
    !************************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Arguments
    integer , intent(in) :: nrow
    integer , intent(in) :: ncol
    real(real64) , intent(in) :: xin(nrow,ncol)
    real(real64) , intent(out) :: xout(ncol,nrow)

    ! Local declarations
    integer :: irow, icol
    !-----------------------------------------------
    !
    xout = transpose(xin)
    !
    return
end subroutine transpose_matrix
