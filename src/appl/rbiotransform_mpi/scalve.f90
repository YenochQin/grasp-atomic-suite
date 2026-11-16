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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
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
