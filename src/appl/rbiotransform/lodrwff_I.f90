      module lodrwff_i
      implicit none (type, external)
      interface
      subroutine lodrwff (name, ntestg)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      integer, intent(in) :: ntestg
      end subroutine
      end interface
      end module
