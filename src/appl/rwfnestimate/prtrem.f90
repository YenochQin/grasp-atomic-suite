

!***********************************************************************
!                                                                      *
      subroutine prtrem(all)
!                                                                      *
!   Prints a list of subshells that remain to be estimated.  ALL  is   *
!   .TRUE.  if all subshells have been estimated and  .FALSE. other-   *
!   wise.                                                              *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use left_c
      use orb_c, only: nw,  np, nak, nh
      use iounit_c
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical , intent(out) :: all
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
!...  /ORB2/
!     COMMON /ORB2/ NCF, NW, PNTRIQ
!     REAL   PNTRIQ
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, iend, ibeg, lenth
      character(len=12) :: cnum
      character(len=80) :: record*80
!-----------------------------------------------
!
!     Print *, ' Entering PRTREM'
!   Determine if there are any subshell radial wavefunctions that
!   remain to be estimated
!
      do i = 1, nw
         if (set(i)) cycle
         all = .false.
         go to 2
      end do
      all = .true.
      go to 4
!
!   Print a list of subshell radial wavefunctions that remain to
!   be estimated; this list is no more than 80 characters wide
!
    2 continue
      write (istde, *) 'The following subshell radial wavefunctions ', &
         'remain to be estimated:'
      iend = 0
      do i = 1, nw
         if (set(i)) cycle
         if (iend > 75) then
            write (istde, *) record(1:iend)
            iend = 0
         endif
         if (iend > 0) then
            ibeg = iend + 1
            iend = ibeg
            record(ibeg:iend) = ' '
         endif
         ibeg = iend + 1
         call convrt (np(i), cnum, lenth)
         iend = ibeg + lenth - 1
         record(ibeg:iend) = cnum(1:lenth)
         if (nak(i) < 0) then
            lenth = 1
         else
            lenth = 2
         endif
         ibeg = iend + 1
         iend = ibeg + lenth - 1
         record(ibeg:iend) = nh(i)(1:lenth)
      end do
      if (iend > 1) write (istde, *) record(1:iend)
!
    4 continue
      return
      end subroutine prtrem
