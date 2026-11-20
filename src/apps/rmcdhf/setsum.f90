!***********************************************************************
      subroutine setsum(file_path)
      use calculation_file_units, only: summary_output_unit
!***********************************************************************
!-----------------------------------------------
!    M O D U L E S
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
      character  :: file_path*(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character(len=9), parameter :: form = "FORMATTED"
      character(len=3), parameter :: status = "NEW"
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: open_status
!-----------------------------------------------

      call openfl (summary_output_unit, file_path, form, status, open_status)
      if (open_status /= 0) then
         write (stderr_unit, *) "Error when opening ", file_path
         error stop
      end if

      return
      end subroutine setsum
