      module solve_i
      implicit none (type, external)
      interface
      subroutine solve (orbital_index, solve_failed, orbital_sign_reversal_count, matching_grid_index, nnp)
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      logical, intent(out) :: solve_failed
      integer, intent(inout) :: orbital_sign_reversal_count
      integer, intent(out) :: matching_grid_index
      integer, intent(out) :: nnp
      end subroutine solve
      end interface
      end module solve_i
