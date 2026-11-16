      module bess_c
      use vast_kind_param, only:  double
      use parameter_def,   only:  max_radial_points
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(double), dimension(2) :: wij
      real(double), dimension(2,2,max_radial_points) :: bessj, bessn
      real(double), dimension(max_radial_points,3) :: bj
      real(double), dimension(max_radial_points) :: tc, td
      end module bess_c
