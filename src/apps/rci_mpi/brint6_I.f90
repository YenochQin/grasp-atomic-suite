      module brint6_i
      implicit none (type, external)
      interface
      subroutine brint6 (orbital_a,orbital_b,orbital_c,orbital_d,radial_integral_order,radial_integral_value)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: orbital_a, orbital_b, orbital_c, orbital_d, radial_integral_order
      real(real64), intent(out) :: radial_integral_value
      end subroutine brint6
      end interface
      end module brint6_i
