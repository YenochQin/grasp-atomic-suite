      module orbord_i
      implicit none (type, external)
      interface
      subroutine orbord (n)
      implicit none (type, external)
      integer :: n
!...Dummy argument N is not referenced in this routine.
      end subroutine orbord
      end interface
      end module orbord_i
