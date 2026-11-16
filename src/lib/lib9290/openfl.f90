!***********************************************************************
!                                                                      *
      subroutine openfl(nfile, filnam, rform, rstat, ierr)
!                                                                      *
!   Issues OPEN for file with unit number NFILE, name FILNAM, format   *
!   RFORM, status RSTAT.  If this is successful the head is position-  *
!   ed to the beginning of the file and IERR is 0; otherwise IERR is   *
!   set to 1.                                                          *
!                                                                      *
!   Call(s) to: [LIB92]:  none.                                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 05 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:04   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use iounit_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nfile
      integer, intent(out) :: ierr
      character (len = *), intent(in) :: filnam
      character (len = *), intent(in) :: rform
      character (len = *), intent(in) :: rstat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ios, loc
!-----------------------------------------------
!
!
      open(nfile, file=filnam, form=rform, status='UNKNOWN', iostat=ios, &
         position='asis')
!
      if (ios == 0) then
         rewind (nfile)
         ierr = 0
      else
         loc = len_trim(filnam)
         write (stderr_unit, *) 'OPENFL: Error opening file ', filnam(1:loc), ' as '&
            , rstat, ';'
         write (stderr_unit, *) 'The argument RSTAT=', rstat, ' is not used !'
         ierr = 1
      endif
!
      return
      end subroutine openfl
