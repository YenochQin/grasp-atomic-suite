      module nuclear_potential_state
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: max_radial_points
      implicit none (type, external)
      private :: real64
      integer :: first_grid_index_beyond_fermi_c
      real(real64) :: negative_radius_times_potential(max_radial_points)
      end module nuclear_potential_state
