      module setcsl_i
      implicit none (type, external)
      interface
      subroutine setcsl (name, ncore, idblk)
      implicit none (type, external)
      character (len = *), intent(in) :: name
      integer :: ncore
      character (len = 8) :: idblk(*)
      end subroutine
      end interface
      end module
