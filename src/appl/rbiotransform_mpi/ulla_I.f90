      module ulla_i
      implicit none (type, external)
      interface
      subroutine ulla (a, u, l, ndim, scr)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ndim
      real(real64), intent(in) :: a(ndim, ndim)
      real(real64), intent(out) :: u(ndim, ndim)
      real(real64), intent(out) :: l(ndim, ndim)
      real(real64), intent(inout) :: scr(*)
      end subroutine
      end interface
      end module
