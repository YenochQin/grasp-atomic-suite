!
!     ------------------------------------------------------------------
!     T R P M A T
!     ------------------------------------------------------------------
!
      subroutine trpmat(xin, nrow, ncol, xout)
!
! XOUT(I,J) = XIN(J,I)
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nrow
      integer , intent(in) :: ncol
      real(real64) , intent(in) :: xin(nrow,ncol)
      real(real64) , intent(out) :: xout(ncol,nrow)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: irow, icol
!-----------------------------------------------
!
      xout = transpose(xin)
!
      return
      end subroutine trpmat
