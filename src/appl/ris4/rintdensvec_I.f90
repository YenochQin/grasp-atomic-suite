      module rintdensvec_i
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine rintdensvec (i, j, dint1vec, nrnuc)
      use vast_kind_param,  only: double
      use parameter_def,    only: nnnw
      use grid_c,           only: n
      real(double), dimension(nnnw,nnnw,n), intent(out) :: dint1vec
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: nrnuc
      end subroutine
      end interface
      end module
