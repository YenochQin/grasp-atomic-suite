      module lulu_i
      implicit none (type, external)
      interface
      subroutine lulu (a, l, u, ndim)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ndim
      real(real64), intent(in) :: a(ndim, ndim)
      real(real64), intent(out) :: l(*)
      real(real64), intent(inout) :: u(*)
      end subroutine lulu
      end interface
      end module lulu_i
