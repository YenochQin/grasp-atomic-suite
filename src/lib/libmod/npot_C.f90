      module npot_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,    only: max_radial_points
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  11:02:52   1/ 2/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: first_grid_index_beyond_fermi_c
      real(real64) :: negative_radius_times_potential(max_radial_points)
      end module npot_c
