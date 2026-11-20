      module in_i
      implicit none (type, external)
      interface
      subroutine in (orbital_index, matching_grid_index, large_component_values, small_component_values, mtp)
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: max_radial_points
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      integer, intent(in) :: matching_grid_index
      real(real64), intent(inout) :: large_component_values(max_radial_points)
      real(real64), intent(inout) :: small_component_values(max_radial_points)
      integer, intent(out) :: mtp
      end subroutine in
      end interface
      end module in_i
