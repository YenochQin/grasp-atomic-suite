!     ------------------------------------------------------------------
!     P R S Y M
!     ------------------------------------------------------------------
!
      subroutine prsym(a, matdim)
! PRINT LOWER HALF OF A SYMMETRIC MATRIX OF DIMENSION MATDIM.
! THE LOWER HALF OF THE MATRIX IS SUPPOSED TO BE IN VECTOR A.
!-----------------------------------------------
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
      integer, intent(in) :: matdim
      real(real64), intent(in) :: a(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: jstart, jstop, i, j
!-----------------------------------------------
      jstart = 1
      jstop = 0
      do i = 1, matdim
         jstart = jstart + i - 1
         jstop = jstop + i
         write (6, 1010) i, (a(j),j=jstart,jstop)
      end do
      return
 1010 format('0',2x,i3,5(1x,e14.7),/,(' ',5x,5(1x,e14.7)))
      return
      end subroutine prsym
