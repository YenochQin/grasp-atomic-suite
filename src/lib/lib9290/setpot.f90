!***********************************************************************
!                                                                      *
      subroutine setpot(j, jp)
!                                                                      *
!   This  subroutine  sets  up the  arrays TF and TG for use by  the   *
!   subprograms IN, OUT, and SBSTEP.                                   *
!                                                                      *
!   The routine uses an algorithm that involves the orbital energy for *
!   finding a join.  When the energy estimate is such that the         *
!   criteria for the join is not met, the program sets the joining     *
!   point in the middle of the grid.                                   *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!      J:  (Input) Index of orbital                                    *
!      JP: (Output) Join point: point where TG changes sign            *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 16 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:38   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use iounit_c
      use def_c
      use grid_c
      use int_c,           only: tf, tg
      use orb_c
      use pote_c,          only: yp
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
      integer, intent(out) :: jp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(double) :: dmhb2c, energy, enefac, rpi, yprpor
      logical :: notset
!-----------------------------------------------
!
!
      notset = .true.
!
!   Define constants
!
      dmhb2c = -h/(2.0d00*c)
      energy = e(j)
      enefac = 2.0d00*c*c - energy
!
!   Set up arrays TF and TG
!
!   Since TF(1) and TG(1) are never used, set them
!   to some arbitrary value
!
      jp = 0
      tf(1) = 0.0d00
      tg(1) = 0.0d00
      do i = 2, n
         rpi = rp(i)
         yprpor = yp(i)*rpor(i)
         tf(i) = dmhb2c*(enefac*rpi + yprpor)
         tg(i) = dmhb2c*(energy*rpi - yprpor)
         if (.not.notset) cycle
         if (abs(sign(1.0d00,tg(i))+sign(1.0d00,tg(1))) >= accy) cycle
         jp = i
         notset = .false.
      end do
!
!   Trap for inappropriate grid
!
      if (jp == 0) then
         write (istde, *) 'SETPOT: Join set to NNNP/2 = ', nnnp/2
         jp = nnnp/2
      endif
!
      return
      end subroutine setpot
