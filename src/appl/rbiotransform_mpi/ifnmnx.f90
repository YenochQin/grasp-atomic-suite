!
!     ------------------------------------------------------------------
!     I F N M N X
!     ------------------------------------------------------------------
!
      integer function ifnmnx (ivec, nel, imxmn)
!
! Smallest or largest value of integer array
!
! IMXMN = 1 => Largest value
! IMXMN = 2 => Smallest
!
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nel
      integer , intent(in) :: imxmn
      integer , intent(in) :: ivec(nel)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ival, iel, ntest
!-----------------------------------------------
!
      ival = ivec(1)
      if (imxmn == 1) then
         ival = max0(maxval(ivec(2:nel)),ival)
      else if (imxmn == 2) then
         ival = min0(minval(ivec(2:nel)),ival)
      else
         write (6, *) ' Stop in IFNMNX '
         write (6, *) ' Improper calue of IMXMN ', imxmn
         stop 'IFNMNX'
      endif
!
      ifnmnx = ival
!
      ntest = 0
      if (ntest /= 0) write (6, *) ' Value returned from IFNMNX', ifnmnx
!
      return
      end function ifnmnx
