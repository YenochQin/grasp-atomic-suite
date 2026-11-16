      module out_i
      implicit none (type, external)
      interface
      subroutine out (j, jp, p, q)
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: max_radial_points
      implicit none (type, external)
      integer, intent(in) :: j
      integer, intent(in) :: jp
      real(real64), intent(inout) :: p(max_radial_points)
      real(real64), intent(inout) :: q(max_radial_points)
      end subroutine out
      end interface
      end module out_i
