      module npot_c
      use vast_kind_param,  only: double
      use parameter_def,    only: max_radial_points
!...Created by Pacific-Sierra Research 77to90  4.3E  11:02:52   1/ 2/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: first_grid_index_beyond_fermi_c
      real(double), dimension(max_radial_points) :: negative_radius_times_potential
      end module npot_c
