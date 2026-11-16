!***********************************************************************
!                                                                      *
      subroutine connect
!                                                                      *
!   The position of an orbital in the merged list is connected to      *
!   the positions in the initial and final state lists                 *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use biorb_c
      use orb_c
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
!-----------------------------------------------
!
!ww      INTEGER PNTRIQ
!     POINTER (PNTRIQ,RIQDUMMY)
!
!
!     Initialize
!
      nnii(:nw) = 0.d0
      nnff(:nw) = 0.d0
!
!   Loop over the orbitals in the merged list
!
      do i = 1, nw
         do j = 1, nwii
            if (np(i)/=npii(j) .or. nak(i)/=nakii(j)) cycle
            nnii(i) = j
         end do
         do j = 1, nwff
            if (np(i)/=npff(j) .or. nak(i)/=nakff(j)) cycle
            nnff(i) = j
         end do
      end do

      return
      end subroutine connect
