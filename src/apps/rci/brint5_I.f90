      module brint5_i
      implicit none (type, external)
      interface
      subroutine brint5 (orbital_a,orbital_b,orbital_c,orbital_d,radial_integral_order,radial_integral_value)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: orbital_a, orbital_b, orbital_c, orbital_d, radial_integral_order
      real(real64), intent(out) :: radial_integral_value
      end subroutine brint5
      end interface
      end module brint5_i
