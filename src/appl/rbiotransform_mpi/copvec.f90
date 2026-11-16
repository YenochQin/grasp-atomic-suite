!
!     ------------------------------------------------------------------
!     C O P V E C
!     ------------------------------------------------------------------
!
      subroutine copvec(from, to, ndim)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ndim
      real(double), dimension(ndim), intent(in)  :: from
      real(double), dimension(ndim), intent(out) :: to
!
      to(:ndim) = from(:ndim)
!
      return
      end subroutine copvec
