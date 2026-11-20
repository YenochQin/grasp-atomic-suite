!
!     ------------------------------------------------------------------
!     C O P V E C
!     ------------------------------------------------------------------
!
      subroutine copvec(from, to, ndim)
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
      real(real64), intent(in)  :: from(ndim)
      real(real64), intent(out) :: to(ndim)
!
      to(:ndim) = from(:ndim)
!
      return
      end subroutine copvec
