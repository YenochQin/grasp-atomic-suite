      module rintdensvec_i
      implicit none (type, external)
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine rintdensvec (first_orbital, second_orbital, dint1vec, nuclear_grid_point_count)
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: max_orbitals
      use radial_grid_state,           only: radial_point_count
      implicit none (type, external)
      real(real64), intent(out) :: dint1vec(max_orbitals, max_orbitals, radial_point_count)
      integer, intent(in) :: first_orbital
      integer, intent(in) :: second_orbital
      integer, intent(in) :: nuclear_grid_point_count
      end subroutine rintdensvec
      end interface
      end module rintdensvec_i
