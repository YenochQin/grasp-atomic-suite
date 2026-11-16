

!***********************************************************************
!                                                                      *
      subroutine setsum
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
!                                                                      *
!   Open the  .sum  file on stream 24.                                 *
!                                                                      *
!   Call(s) to: [LIB92]:  OPENFL.                                      *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i, only: openfl
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ierr
      character :: filnam*256, defnam*11, form*11, status*3
!-----------------------------------------------
!
!   File  erwf.sum  is FORMATTED
!
      defnam = "erwf.sum"
      form = "FORMATTED"
      status = "NEW"
!
      write (stderr_unit, *) "File  erwf.sum  will be created as the ", &
         "ERWF SUMmary File; "
      write (stderr_unit, *) "enter another file name if this is not ", &
         "acceptable; null otherwise:"
      read (*, "(A)") filnam
!
      if (len_trim(filnam) == 0) filnam = defnam
!
    1 continue
      call openfl (24, filnam, form, status, ierr)
      if (ierr /= 0) then
    2    continue
         write (stderr_unit, *) "Enter a name for the ERWF SUMmary File ", &
            "that is to be created:"
         read (*, "(A)") filnam
         if (len_trim(filnam) == 0) goto 2
         goto 1
      end if
!
      return
      end subroutine setsum
