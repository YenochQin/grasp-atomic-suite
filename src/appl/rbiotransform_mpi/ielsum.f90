!
!     ------------------------------------------------------------------
!     I E L S U M
!     ------------------------------------------------------------------
!
      integer function ielsum (ivec, nelmnt)
!
! Sum elements of integer array
!
!************************************************************************
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
      integer, intent(in) :: nelmnt
      integer, intent(in) :: ivec(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: isum, iel
!-----------------------------------------------
!
      isum = 0
      isum = sum(ivec(:nelmnt))
!
      ielsum = isum
!
      return
      end function ielsum
