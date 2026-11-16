!***********************************************************************
!                                                                      *
      subroutine setcsl
!                                                                      *
!   Open, check, load data from and close the  .csl  file. This file   *
!   is always attached to stream 21.                                   *
!                                                                      *
!   Call(s) to: [RCI92]: LENGTH, LODCSL, OPENFL.                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      use lodcsl_i
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ierr, ios, ncore
      logical :: found
      character :: filnam*256, record*15, defnam*11, form*11, status*3
!-----------------------------------------------
!
!
!   The  .csl  file is FORMATTED; it must exist
!
      defnam = 'rcsl.inp'
      form = 'FORMATTED'
      status = 'OLD'
!
!   Look for  grasp92.csl
!
      inquire(file=defnam, exist=found)
!
      if (found) then
         filnam = defnam
      else
         write (6, *) 'rcsl.inp does not exist'
         stop
      endif
!
      call openfl (21, filnam, form, status, ierr)
      if (ierr == 1) then
         write (6, *) 'Error when opening rcsl.inp'
         stop
      endif
!
!   Check the first record of the file; if not as expected, try again
!
      read (21, '(1A15)', iostat=ios) record
      if (ios/=0 .or. record(1:15)/='Core subshells:') then
         write (6, *) 'Not a Configuration Symmetry List File;'
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
      end subroutine setcsl
