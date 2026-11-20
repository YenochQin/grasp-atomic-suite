      module setmix_i
      implicit none (type, external)
      interface
      subroutine setmix (name)
      implicit none (type, external)
      character (len = *), intent(in) :: name
      end subroutine setmix
      end interface
      end module setmix_i
