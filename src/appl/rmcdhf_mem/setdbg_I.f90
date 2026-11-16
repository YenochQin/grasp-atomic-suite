      module setdbg_i
      implicit none (type, external)
      interface
      subroutine setdbg (dbgfile)
      implicit none (type, external)
      character (len = *), intent(in) :: dbgfile
      end subroutine
      end interface
      end module
