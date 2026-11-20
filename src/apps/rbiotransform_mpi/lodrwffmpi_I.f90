      module lodrwff_i
      implicit none (type, external)
      interface
      subroutine lodrwff (state_name, ntestg)
      implicit none (type, external)
      character (len = 128), intent(in) :: state_name
      integer, intent(in) :: ntestg
      end subroutine lodrwff
      end interface
      end module lodrwff_i
