!
!     ------------------------------------------------------------------
!     I N P R O D
!     ------------------------------------------------------------------
!
      real(kind(0.0d0)) function inprod (a, b, ndim)
!      CALCULATE SCALAR PRODUCT BETWEEN TO VECTORS A,B
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
      real(double), intent(in) :: a(*)
      real(double), intent(in) :: b(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
!
      inprod = dot_product(a(:ndim),b(:ndim))
!
      return
      end function inprod
