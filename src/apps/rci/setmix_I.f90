      module setmix_i
      implicit none (type, external)
      interface
      subroutine setmix (name, block_labels)
      implicit none (type, external)
      character (len = *), intent(in) :: name
      character (len = 8) :: block_labels(*)
      end subroutine setmix
      end interface
      end module setmix_i
