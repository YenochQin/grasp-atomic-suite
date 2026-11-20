      module mohr_i
      implicit none (type, external)
      interface
      subroutine mohr (principal_n, kappa, nuclear_charge, self_energy_factor)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, parameter :: numval = 12
      integer, intent(in) :: principal_n
      integer, intent(in) :: kappa
      real(real64) :: nuclear_charge
      real(real64), intent(out) :: self_energy_factor
      end subroutine mohr
      end interface
      end module mohr_i
