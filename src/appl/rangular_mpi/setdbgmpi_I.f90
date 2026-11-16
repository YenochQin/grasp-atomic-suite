      module setdbgmpi_i
      implicit none (type, external)
      interface
      subroutine setdbgmpi (debug, fullname)
      implicit none (type, external)
      logical :: debug
      character(len=*) :: fullname
      end subroutine setdbgmpi
      end interface
      end module setdbgmpi_i
