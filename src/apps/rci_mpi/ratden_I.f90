      module ratden_i
      implicit none (type, external)
      interface
      real(real64) function ratden (large_component_values, small_component_values, last_radial_index, principal_n, kappa, &
         & effective_nuclear_charge)
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_radial_points
      implicit none (type, external)
      real(real64), intent(in) :: large_component_values(max_radial_points)
      real(real64), intent(in) :: small_component_values(max_radial_points)
      integer, intent(in) :: last_radial_index
      integer :: principal_n
      integer :: kappa
      real(real64) :: effective_nuclear_charge
      end function ratden
      end interface
      end module ratden_i
