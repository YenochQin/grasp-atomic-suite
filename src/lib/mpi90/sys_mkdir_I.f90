      module sys_mkdir_i
      implicit none (type, external)
      interface
      subroutine sys_mkdir(dir, lendir, ierr)
      implicit none (type, external)
      character(len=*), intent(in):: dir
      integer, intent(in):: lendir
      integer, intent(out):: ierr
      end subroutine
      end interface
      end module
