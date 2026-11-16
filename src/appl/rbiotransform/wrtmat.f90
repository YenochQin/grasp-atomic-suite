!***********************************************************************
!                                                                      *
      subroutine wrtmat(a, nrow, ncol, nmrow, nmcol)
!                                                                      *
!***********************************************************************
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
      integer, intent(in) :: nrow
      integer, intent(in) :: ncol
      integer, intent(in) :: nmrow
      integer, intent(in) :: nmcol
      real(real64), intent(in) :: a(nmrow,nmcol)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
!-----------------------------------------------

      do i = 1, nrow
         write (6, 1010) i, (a(i,j),j=1,ncol)
      end do

 1010 format('0',i5,2x,4(1x,e14.7),/,(' ',7x,4(1x,e14.7)))

      return
      end subroutine wrtmat
