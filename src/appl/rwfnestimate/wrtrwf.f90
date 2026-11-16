

!***********************************************************************
!                                                                      *
      subroutine wrtrwf
!                                                                      *
!   Open, write a header and all subshell radial wavefunctions, and    *
!   close the  .rwf  file.                                             *
!                                                                      *
!   Call(s) to: [LIB92]: DALLOC, OPENFL.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use memory_man
      use grid_c
      use orb_c
      use wave_c, only: pz, pf, qf, mf
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: newunit, ierr, j, mfj, i
      logical :: imopened
      character :: filnam*128
!-----------------------------------------------
!
!

      filnam = 'rwfn.inp'

      do newunit = 23, 99                        ! 23 is a historical value
         inquire(unit=newunit, opened=imopened)
         if (imopened) cycle
         exit                                    ! should be the normal exit point
      end do

      if (newunit == 100) then
         write (istde, *) 'All unit numbers from 23 to 99 are BUSY!'
         stop
      endif

      call openfl (newunit, filnam, 'UNFORMATTED', 'NEW', ierr)
      if (ierr == 1) then
         write (istde, *) 'Error when opening "', filnam(1:len_trim(filnam)), &
            '"'

         stop
      endif
!
!   Write the file header
!
      write (newunit) 'G92RWF'
!
!   Write out the radial wavefunctions
!
      do j = 1, nw
         mfj = mf(j)
         write (newunit) np(j), nak(j), e(j), mfj
         write (newunit) pz(j), (pf(i,j),i=1,mfj), (qf(i,j),i=1,mfj)
         write (newunit) (r(i),i=1,mfj)
      end do

      close(newunit)
!
!   Deallocate the storage for the radial wavefunctions
!
      call dalloc (pf, 'PF', 'WRTRWF')
      call dalloc (qf, 'QF', 'WQRTRWF')
!
      return
      end subroutine wrtrwf
