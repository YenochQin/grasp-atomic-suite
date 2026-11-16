      module lodrwfi_i
      implicit none (type, external)
      interface
      subroutine lodrwfi (name, ntestg)
      implicit none (type, external)
      character (len = 128), intent(in) :: name
      integer, intent(in) :: ntestg
      end subroutine
      end interface
      end module
