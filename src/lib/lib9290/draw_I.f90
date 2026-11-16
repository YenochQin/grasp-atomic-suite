      module draw_i
      implicit none (type, external)
      interface
      subroutine draw (p, sp, q, sq, mf)
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def, only: max_radial_points
      implicit none (type, external)
      real(real64), intent(in) :: p(max_radial_points)
      real(real64), intent(in) :: sp
      real(real64), intent(in) :: q(max_radial_points)
      real(real64), intent(in) :: sq
      integer, intent(in) :: mf
!VAST.../GRID/ R(IN), N(IN)
!...This routine performs I/O.
      end subroutine
      end interface
      end module
