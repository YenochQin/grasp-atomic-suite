      module klamaq_i
      implicit none (type, external)
      interface
      subroutine klamaq (principal_n, kappa, nuclear_charge, self_energy_factor)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: principal_n
      integer, intent(in) :: kappa
      real(real64), intent(in) :: nuclear_charge
      real(real64), intent(out) :: self_energy_factor
      end subroutine klamaq
      end interface
      end module klamaq_i
