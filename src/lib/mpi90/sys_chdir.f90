module sys_chdir_i
    implicit none (type, external)
    private
    public :: sys_chdir
contains
!***********************************************************************
      subroutine sys_chdir (dir, lendir, ierr)
!
!   This routine changes current working directory to dir.
!   lendir is the length of character string dir;
!   ierr will be zero if successful, otherwise non-zero;
!
!   Xinghong He  98-08-21
!
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use suite_system, only: change_directory
      implicit none (type, external)
      character(len=*), intent(in):: dir
      integer, intent(in):: lendir
      integer, intent(out):: ierr


      ierr = change_directory (dir(1:lendir))

      return
      end subroutine sys_chdir
end module sys_chdir_i
