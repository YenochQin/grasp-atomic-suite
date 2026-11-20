      module setsum_i
      implicit none (type, external)
      interface
      subroutine setsum (name)
      implicit none (type, external)
      character (len = *), intent(in) :: name
      end subroutine setsum
      end interface
      end module setsum_i
