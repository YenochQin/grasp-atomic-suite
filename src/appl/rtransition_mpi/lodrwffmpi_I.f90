      module lodrwff_i
      implicit none (type, external)
      interface
      subroutine lodrwff (name)
      implicit none (type, external)
      character (len = 128), intent(in) :: name
      end subroutine lodrwff
      end interface
      end module lodrwff_i
