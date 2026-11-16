      module setmix_i
      implicit none (type, external)
      interface
      subroutine setmix (name, idblk)
      implicit none (type, external)
      character (len = *), intent(in) :: name
      character (len = 8) :: idblk(*)
      end subroutine setmix
      end interface
      end module setmix_i
