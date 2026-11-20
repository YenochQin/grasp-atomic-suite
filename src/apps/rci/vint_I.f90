      module vint_i
      implicit none (type, external)
      interface
      subroutine vint (first_orbital,second_orbital,radial_integral_value)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: first_orbital, second_orbital
      real(real64), intent(out) :: radial_integral_value
      end subroutine vint
      end interface
      end module vint_i
