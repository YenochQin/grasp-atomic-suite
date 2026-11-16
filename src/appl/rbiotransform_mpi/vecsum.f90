!
!     ------------------------------------------------------------------
!     V E C S U M
!     ------------------------------------------------------------------
!
      subroutine vecsum(c, a, b, faca, facb, ndim)
!
!     CACLULATE THE VECTOR C(I)=FACA*A(I)+FACB*B(I)
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ndim
      real(real64), intent(in) :: faca
      real(real64), intent(in) :: facb
      real(real64), intent(in)  :: a(ndim)
      real(real64), intent(in)  :: b(ndim)
      real(real64), intent(out) :: c(ndim)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(real64) :: s
!-----------------------------------------------
!
      if (faca/=0.0_real64 .and. facb/=0.0_real64) then
         c(:ndim) = faca*a(:ndim) + facb*b(:ndim)
!
      else if (faca==0.0_real64 .and. facb/=0.0_real64) then
         c(:ndim) = facb*b(:ndim)
!
      else if (faca/=0.0_real64 .and. facb==0.0_real64) then
         c(:ndim) = faca*a(:ndim)
!
      else if (faca==0.0_real64 .and. facb==0.0_real64) then
         c(:ndim) = 0.0_real64
      endif
!
      return
      end subroutine vecsum
