!***********************************************************************
!                                                                      *
      subroutine setcslb(name, ncore, igg)
!-----------------------------------------------
!                                                                      *
!   Open, check, load data from and close the  .csl  file. This file   *
!   is always attached to stream 21.                                   *
!                                                                      *
!   Call(s) to: [RCI92]: LENGTH, LODCSL, OPENFL.                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      use lodcslbio_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ncore, igg
      character, intent(in) :: name*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, ierr, ios
      logical :: found
      character :: filnam*256, record*15, defnam*11, form*11, status*3
!
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
         write (6, *) 'Error when opening', filnam
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
      call lodcslbio (ncore,igg)
!
!   Close the  .csl  file
!
      close(21)
!
      return
      end subroutine setcslb
