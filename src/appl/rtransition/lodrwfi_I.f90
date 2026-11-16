      module lodrwfi_i
      implicit none (type, external)
      interface
      subroutine lodrwfi (name)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      end subroutine lodrwfi
      end interface
      end module lodrwfi_i
