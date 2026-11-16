!
!     ------------------------------------------------------------------
!     S E T V E C
!     ------------------------------------------------------------------
!
      subroutine setvec(vector, value, ndim)
!
! VECTOR (*) = VALUE
!
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
      real(real64), intent(in) :: value
      real(real64), intent(out) :: vector(ndim)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
!
      vector(:ndim) = value
!
      return
      end subroutine setvec
