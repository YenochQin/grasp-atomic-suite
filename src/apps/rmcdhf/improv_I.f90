      module improv_i
      implicit none (type, external)
      interface
      subroutine improv (extended_optimal_level, orbital_index, lsort, dampmx)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      logical, intent(in) :: extended_optimal_level
      integer, intent(in) :: orbital_index
      logical :: lsort
      real(real64), intent(inout) :: dampmx
      end subroutine improv
      end interface
      end module improv_i
