      module setcsl_i
      implicit none (type, external)
      interface
      subroutine setcsl (name, ncore, nblkin, idblk)
      implicit none (type, external)
      character (len = *) :: name
      integer :: ncore
      integer, intent(in) :: nblkin
      character (len = 8) :: idblk(*)
      end subroutine setcsl
      end interface
      end module setcsl_i
