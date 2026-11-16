      module setsum_i
      implicit none (type, external)
      interface
      subroutine setsum (name)
      implicit none (type, external)
      character (len = *), intent(in) :: name
      end subroutine
      end interface
      end module
