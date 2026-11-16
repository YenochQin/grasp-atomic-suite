      module lodrwfi_i
      implicit none (type, external)
      interface
      subroutine lodrwfi (name)
      implicit none (type, external)
      character (len = 128), intent(in) :: name
      end subroutine
      end interface
      end module
