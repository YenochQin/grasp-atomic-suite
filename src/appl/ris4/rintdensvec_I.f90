      module rintdensvec_i
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine rintdensvec (i, j, dint1vec, nrnuc)
      use vast_kind_param,  only: double
      use parameter_def,    only: max_orbitals
      use grid_c,           only: radial_point_count
      real(double), dimension(max_orbitals,max_orbitals,radial_point_count), intent(out) :: dint1vec
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: nrnuc
      end subroutine
      end interface
      end module
