      module setdbgmpi_i
      implicit none (type, external)
      interface
      subroutine setdbgmpi (debug_file_path)
      implicit none (type, external)
      character (len = *), intent(in) :: debug_file_path
      end subroutine setdbgmpi
      end interface
      end module setdbgmpi_i
