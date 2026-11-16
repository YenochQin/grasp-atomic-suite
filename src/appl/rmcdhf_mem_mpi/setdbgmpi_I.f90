      module setdbgmpi_i
      implicit none (type, external)
      interface
      subroutine setdbgmpi (dbgfile)
      implicit none (type, external)
      character (len = *), intent(in) :: dbgfile
      end subroutine
      end interface
      end module
