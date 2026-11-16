!
!     ------------------------------------------------------------------
!     V E C S U M
!     ------------------------------------------------------------------
!
      subroutine vecsum(c, a, b, faca, facb, ndim)
!
!     CACLULATE THE VECTOR C(I)=FACA*A(I)+FACB*B(I)
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
      integer, intent(in) :: ndim
      real(double), intent(in) :: faca
      real(double), intent(in) :: facb
      real(double), dimension(ndim), intent(in)  :: a
      real(double), dimension(ndim), intent(in)  :: b
      real(double), dimension(ndim), intent(out) :: c
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(double) :: s
!-----------------------------------------------
!
      if (faca/=0.0d0 .and. facb/=0.0d0) then
         c(:ndim) = faca*a(:ndim) + facb*b(:ndim)
!
      else if (faca==0.0d0 .and. facb/=0.0d0) then
         c(:ndim) = facb*b(:ndim)
!
      else if (faca/=0.0d0 .and. facb==0.0d0) then
         c(:ndim) = faca*a(:ndim)
!
      else if (faca==0.0d0 .and. facb==0.0d0) then
         c(:ndim) = 0.0d0
      endif
!
      return
      end subroutine vecsum
