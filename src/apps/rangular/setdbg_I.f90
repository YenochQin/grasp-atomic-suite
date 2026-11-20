      module setdbg_i
      implicit none (type, external)
      interface
      subroutine setdbg (debug, file_path)
      implicit none (type, external)
      logical :: debug
      character(len=*) :: file_path
      end subroutine setdbg
      end interface
      end module setdbg_i
