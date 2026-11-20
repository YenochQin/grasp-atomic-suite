      module wrtrwf_i
      implicit none (type, external)
      interface
      subroutine wrtrwf(file_path)
      implicit none (type, external)
      character(len=*), intent(in) :: file_path
      end subroutine wrtrwf
      end interface
      end module wrtrwf_i
