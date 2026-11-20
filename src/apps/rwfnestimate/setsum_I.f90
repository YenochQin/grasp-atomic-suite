      module setsum_i
      implicit none (type, external)
      interface
      subroutine setsum(default_file_path)
      implicit none (type, external)
      character(len=*), intent(in) :: default_file_path
      end subroutine setsum
      end interface
      end module setsum_i
