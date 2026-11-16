      module orbord_i
      implicit none (type, external)
      interface
      subroutine orbord (n)
      implicit none (type, external)
      integer :: n
!...Dummy argument N is not referenced in this routine.
      end subroutine
      end interface
      end module
