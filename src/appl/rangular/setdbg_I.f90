      module setdbg_i
      implicit none (type, external)
      interface
      subroutine setdbg (debug, fullname)
      implicit none (type, external)
      logical :: debug
      character(len=*) :: fullname
      end subroutine setdbg
      end interface
      end module
