      module lodrwff_i
      implicit none (type, external)
      interface
      subroutine lodrwff (state_name)
      implicit none (type, external)
      character (len = 128), intent(in) :: state_name
      end subroutine lodrwff
      end interface
      end module lodrwff_i
