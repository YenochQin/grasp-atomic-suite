      module lodrwfi_i
      implicit none (type, external)
      interface
      subroutine lodrwfi (state_name, ntestg)
      implicit none (type, external)
      character (len = 24), intent(in) :: state_name
      integer, intent(in) :: ntestg
      end subroutine lodrwfi
      end interface
      end module lodrwfi_i
