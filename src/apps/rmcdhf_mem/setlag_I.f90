      module setlag_i
      implicit none (type, external)
      interface
      subroutine setlag (extended_optimal_level)
      implicit none (type, external)
      logical :: extended_optimal_level
!...This routine performs I/O.
      end subroutine setlag
      end interface
      end module setlag_i
