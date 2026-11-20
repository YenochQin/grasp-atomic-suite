      module wrtrwf_i
      implicit none (type, external)
      interface
      subroutine wrtrwf(filnam)
      implicit none (type, external)
      character(len=*), intent(in) :: filnam
      end subroutine wrtrwf
      end interface
      end module wrtrwf_i
