      module setsum_i
      implicit none (type, external)
      interface
      subroutine setsum (file_path)
      implicit none (type, external)
      character (len = *), intent(in) :: file_path
      end subroutine setsum
      end interface
      end module setsum_i
