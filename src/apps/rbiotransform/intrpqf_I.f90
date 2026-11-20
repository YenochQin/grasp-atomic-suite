      module intrpqf_i
      implicit none (type, external)
      interface
      subroutine intrpqf (source_large_component, source_small_component, source_grid_size, source_radial_grid, orbital_index, &
         & original_norm_squared)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: source_large_component(*)
      real(real64), intent(in) :: source_small_component(*)
      integer, intent(in) :: source_grid_size
      real(real64), intent(in) :: source_radial_grid(*)
      integer, intent(in) :: orbital_index
      real(real64), intent(out) :: original_norm_squared
      end subroutine intrpqf
      end interface
      end module intrpqf_i
