      module newe_i
      implicit none (type, external)
      interface
      subroutine newe (orbital_index, first_oscillation_sign, effective_principal_n, node_count_difference, energy_correction, &
         & solve_failed, orbital_sign_reversal_count)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      real(real64), intent(in) :: first_oscillation_sign
      integer, intent(in) :: effective_principal_n
      integer, intent(in) :: node_count_difference
      real(real64), intent(inout) :: energy_correction
      logical, intent(out) :: solve_failed
      integer, intent(out) :: orbital_sign_reversal_count
      end subroutine newe
      end interface
      end module newe_i
