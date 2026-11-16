!
!
!     ------------------------------------------------------------------
!     S C A L V E
!     ------------------------------------------------------------------
!
      subroutine scalve(vector, factor, ndim)
!
! CALCULATE SCALAR(FACTOR) TIMES VECTOR
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
      real(real64), intent(in) :: factor
      real(real64), intent(inout) :: vector(ndim)
!-----------------------------------------------
!
      vector(:ndim) = vector(:ndim)*factor
!
      return
      end subroutine scalve
