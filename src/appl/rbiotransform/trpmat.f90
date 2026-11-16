!
!     ------------------------------------------------------------------
!     T R P M A T
!     ------------------------------------------------------------------
!
      subroutine trpmat(xin, nrow, ncol, xout)
!
! XOUT(I,J) = XIN(J,I)
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nrow
      integer , intent(in) :: ncol
      real(double) , intent(in) :: xin(nrow,ncol)
      real(double) , intent(out) :: xout(ncol,nrow)
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
