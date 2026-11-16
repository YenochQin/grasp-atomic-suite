      module sys_chdir_i
      implicit none (type, external)
      interface
      subroutine sys_chdir (dir, lendir, ierr)
      implicit none (type, external)
      character(len=*), intent(in):: dir
      integer, intent(in):: lendir
      integer, intent(out):: ierr
      end subroutine
      end interface
      end module
