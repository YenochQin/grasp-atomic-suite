!***********************************************************************
!                                                                      *
      subroutine tcsl(n)
!                                                                      *
!   This subroutine transfers data to the initial and final state      *
!   common blocks                                                      *
!                                                                      *
!   Written by Per Jonsson                Last revision: June   1996   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
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
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
!
!  Initial state commons
!
!
!  Final state commons
!
!
      if (n == 1) then
         ncfii = ncf
         nwii = nw
         npii(:nw) = np(:nw)
         nakii(:nw) = nak(:nw)
         nklii(:nw) = nkl(:nw)
         nkjii(:nw) = nkj(:nw)
         nhii(:nw) = nh(:nw)
      else
         ncfff = ncf
         nwff = nw
         npff(:nw) = np(:nw)
         nakff(:nw) = nak(:nw)
         nklff(:nw) = nkl(:nw)
         nkjff(:nw) = nkj(:nw)
         nhff(:nw) = nh(:nw)
      endif

      return
      end subroutine tcsl
