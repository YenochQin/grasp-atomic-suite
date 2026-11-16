!***********************************************************************
!                                                                      *
      subroutine setcsla(name, ncore)
!                                                                      *
!   Open, check, load data from and close the  .csl  file. This file   *
!   is always attached to stream 21.                                   *
!                                                                      *
!   Call(s) to: [RCI92]: LODCSL, OPENFL.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:33   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      use lodcsl_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ncore
      character (len = 24), intent(in) :: name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer   :: k, ierr, ios
      character(len=3)   :: status
      character(len=9)   :: form
      character(len=15)  :: record
      character(len=256) :: filnam
!-----------------------------------------------
!
!   The  .csl  file is FORMATTED; it must exist
!
      k = index(name,' ')
      filnam = name(1:k-1)//'.c'
      form = 'FORMATTED'
      status = 'OLD'

!
      call openfl (21, filnam, form, status, ierr)
      if (ierr == 1) then
         write (stderr_unit, *) 'Error when opening', filnam
         stop
      endif
!
!   Check the first record of the file; if not as expected, try again
!
      read (21, '(1A15)', iostat=ios) record
      if (ios/=0 .or. record(1:15)/='Core subshells:') then
         write (stderr_unit, *) 'Not a Configuration Symmetry List File;'
         close(21)
      endif
!
!   Load data from the  .csl  file
!
      call lodcsl (ncore)
!
!   Close the  .csl  file
!
      close(21)
!
      return
      end subroutine setcsla
