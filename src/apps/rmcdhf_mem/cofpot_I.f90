      module cofpot_i
      implicit none (type, external)
      interface
      subroutine cofpot (extended_optimal_level, orbital_index, radial_point_count)
      implicit none (type, external)
      logical :: extended_optimal_level
      integer :: orbital_index
      integer :: radial_point_count
      end subroutine cofpot
      end interface
      end module cofpot_i
