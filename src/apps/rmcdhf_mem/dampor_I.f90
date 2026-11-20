      module dampor_i
      implicit none (type, external)
      interface
      subroutine dampor (orbital_index, orbital_sign_reversal_count, orbital_damping_factor)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      integer, intent(inout) :: orbital_sign_reversal_count
      real(real64), intent(in) :: orbital_damping_factor
      end subroutine dampor
      end interface
      end module dampor_i
