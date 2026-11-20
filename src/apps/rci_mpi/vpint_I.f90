      module vpint_i
      implicit none (type, external)
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine vpint (first_orbital,second_orbital,radial_integral_value)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer , intent(inout) :: first_orbital, second_orbital
      real(real64), intent(out) :: radial_integral_value
      end subroutine vpint
      end interface
      end module vpint_i
