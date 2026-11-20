      module setdbg_i
      implicit none (type, external)
      interface
      subroutine setdbg(default_file_path)
      implicit none (type, external)
      character(len=*), intent(in) :: default_file_path
      end subroutine setdbg
      end interface
      end module setdbg_i
