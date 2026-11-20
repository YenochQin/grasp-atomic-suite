!***********************************************************************
!                                                                      *
      subroutine setsum(file_path)
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use file_open, only: openfl
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=*), intent(in) :: file_path
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: open_status
      character(len=120) :: filnam  ! The original did not compile
      character(len=11) :: form
      character(len=3)  :: status
!-----------------------------------------------

      form = "FORMATTED"
      status = "NEW"
!
      write (stderr_unit, *) "File ", file_path, " will be created as the", &
         " GENMCP SUMmary File;"
      write (stderr_unit, *) "enter another file name if this is not ", &
         "acceptable; null otherwise:"
      read (*, "(A)") filnam
!
      if (len_trim(filnam) == 0) filnam = file_path
!
    1 continue
      call openfl (24, filnam, form, status, open_status)
      if (open_status /= 0) then
    2    continue
         write (stderr_unit, *) "Enter a name for the GENMCP SUMmary", &
            " File that is to be created:"
         read (*, "(A)") filnam
         if (len_trim(filnam) == 0) goto 2
         goto 1
      end if
!
      return
      end subroutine setsum
