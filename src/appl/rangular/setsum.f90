!***********************************************************************
!                                                                      *
      subroutine setsum(fullname)
!                                                                      *
!   Open the  .sum  file on stream 24.                                 *
!                                                                      *
!   Call(s) to: [LIB92]:  OPENFL.                                      *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 11 Nov 1992   *
!   Modified by Xinghong He               Last revision:  3 Jul 1998   *
!
!  File shared by mcpblk, mcpmpi
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:54:22   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=*), intent(in) :: fullname
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ierr
      character(len=120) :: filnam  ! The original did not compile
      character(len=11) :: form
      character(len=3)  :: status
!-----------------------------------------------

      form = 'FORMATTED'
      status = 'NEW'
!
      write (istde, *) 'File ', fullname, ' will be created as the', &
         ' GENMCP SUMmary File;'
      write (istde, *) 'enter another file name if this is not ', &
         'acceptable; null otherwise:'
      read (*, '(A)') filnam
!
      if (len_trim(filnam) == 0) filnam = fullname
!
    1 continue
      call openfl (24, filnam, form, status, ierr)
      if (ierr /= 0) then
    2    continue
         write (istde, *) 'Enter a name for the GENMCP SUMmary', &
            ' File that is to be created:'
         read (*, '(A)') filnam
         if (len_trim(filnam) == 0) go to 2
         go to 1
      endif
!
      return
      end subroutine setsum
