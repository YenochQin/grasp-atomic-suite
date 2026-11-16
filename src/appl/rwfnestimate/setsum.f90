

!***********************************************************************
!                                                                      *
      subroutine setsum
      use iounit_c
!                                                                      *
!   Open the  .sum  file on stream 24.                                 *
!                                                                      *
!   Call(s) to: [LIB92]:  OPENFL.                                      *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ierr
      character :: filnam*256, defnam*11, form*11, status*3
!-----------------------------------------------
!
!   File  erwf.sum  is FORMATTED
!
      defnam = 'erwf.sum'
      form = 'FORMATTED'
      status = 'NEW'
!
      write (stderr_unit, *) 'File  erwf.sum  will be created as the ', &
         'ERWF SUMmary File; '
      write (stderr_unit, *) 'enter another file name if this is not ', &
         'acceptable; null otherwise:'
      read (*, '(A)') filnam
!
      if (len_trim(filnam) == 0) filnam = defnam
!
    1 continue
      call openfl (24, filnam, form, status, ierr)
      if (ierr /= 0) then
    2    continue
         write (stderr_unit, *) 'Enter a name for the ERWF SUMmary File ', &
            'that is to be created:'
         read (*, '(A)') filnam
         if (len_trim(filnam) == 0) go to 2
         go to 1
      endif
!
      return
      end subroutine setsum
