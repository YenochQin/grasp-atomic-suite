      module setcof_i
      implicit none (type, external)
      interface
      subroutine setcof (extended_optimal_level, orbital_index)
      implicit none (type, external)
      logical :: extended_optimal_level
      integer, intent(in) :: orbital_index
      end subroutine setcof
      end interface
      end module setcof_i
