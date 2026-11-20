      module setlagmpi_i
      implicit none (type, external)
      interface
      subroutine setlagmpi (extended_optimal_level)
      implicit none (type, external)
      logical :: extended_optimal_level
!...This routine performs I/O.
      end subroutine setlagmpi
      end interface
      end module setlagmpi_i
