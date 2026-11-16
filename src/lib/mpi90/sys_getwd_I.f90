      module sys_getwd_i
      implicit none (type, external)
      interface
      subroutine sys_getwd(dir,lcwd)
      implicit none (type, external)
      character(len=*), intent(out):: dir
      integer, intent(out):: lcwd
      end subroutine
      end interface
      end module
