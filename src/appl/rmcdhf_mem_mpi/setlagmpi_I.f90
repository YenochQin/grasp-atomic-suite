      module setlagmpi_i
      implicit none (type, external)
      interface
      subroutine setlagmpi (eol)
      implicit none (type, external)
      logical :: eol
!...This routine performs I/O.
      end subroutine
      end interface
      end module
