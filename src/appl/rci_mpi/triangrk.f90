!***********************************************************************
!                                                                      *
      logical function triangrk (la, k, lb)
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: la
      integer, intent(in) :: k
      integer, intent(in) :: lb
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!   Perform the triangularity check
!
      if (mod(k + la + lb,2) /= 0) then
         triangrk = .false.
      else
         if (abs(la - lb) > k) then
            triangrk = .false.
         else if (la + lb < k) then
            triangrk = .false.
         else
            triangrk = .true.
         endif
      endif

      return
      end function triangrk
