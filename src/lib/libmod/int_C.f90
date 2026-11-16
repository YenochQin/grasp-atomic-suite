      module int_c
      use vast_kind_param,  only: double
      use parameter_def,    only: max_radial_points
!...Created by Pacific-Sierra Research 77to90  4.3E  06:37:37  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: mtp0
      real(double), dimension(max_radial_points) :: p, q
      real(double) :: p0, q0
      real(double), dimension(max_radial_points) :: tf, tg, xu, xv
      end module int_c
