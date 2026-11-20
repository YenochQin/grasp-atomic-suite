      module brint2_i
      implicit none (type, external)
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine brint2 (orbital_a,orbital_b,orbital_c,orbital_d,k,radial_integral_value)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: orbital_a, orbital_b, orbital_c, orbital_d, k
      real(real64), intent(out) :: radial_integral_value
      end subroutine brint2
      end interface
      end module brint2_i
