!***********************************************************************
!                                                                      *
      subroutine maxarr(j)
!                                                                      *
!   This subroutine finds the least self-consistent orbital            *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 10 Dec 1992   *
!                                                                      *
! J initialized to zero
! XHH 1997.02.14
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use fixd_c
      use orb_c
      use scf_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(out) :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(double) :: dlrgst
!-----------------------------------------------
!
!
      j = 0
      dlrgst = 0.d0
      do i = 1, nw
         if (lfix(i)) cycle
         if (scnsty(i) <= dlrgst) cycle
         dlrgst = scnsty(i)
         j = i
      end do
!
      return
      end subroutine maxarr
