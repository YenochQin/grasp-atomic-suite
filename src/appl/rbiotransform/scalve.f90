!
!
!     ------------------------------------------------------------------
!     S C A L V E
!     ------------------------------------------------------------------
!
      subroutine scalve(vector, factor, ndim)
!
! CALCULATE SCALAR(FACTOR) TIMES VECTOR
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
      real(double), intent(in) :: factor
      real(double), dimension(ndim), intent(inout) :: vector
!-----------------------------------------------
!
      vector(:ndim) = vector(:ndim)*factor
!
      return
      end subroutine scalve
