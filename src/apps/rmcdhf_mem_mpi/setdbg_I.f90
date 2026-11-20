      module setdbg_i
      implicit none (type, external)
      interface
      subroutine setdbg (debug_file_path)
      implicit none (type, external)
      character (len = *), intent(in) :: debug_file_path
      end subroutine setdbg
      end interface
      end module setdbg_i
