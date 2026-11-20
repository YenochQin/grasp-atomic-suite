      module iabint_i
      implicit none (type, external)
      interface
      subroutine iabint (first_orbital,second_orbital,radial_integral_value)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(inout) :: first_orbital, second_orbital
      real(real64), intent(out) :: radial_integral_value
      end subroutine iabint
      end interface
      end module iabint_i
