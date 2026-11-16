      module setcsl_i
      implicit none (type, external)
      interface
      subroutine setcsl (name, ncore, idblk)
      implicit none (type, external)
      character (len = *), intent(in) :: name
      integer :: ncore
      character (len = 8) :: idblk(*)
      end subroutine setcsl
      end interface
      end module setcsl_i
