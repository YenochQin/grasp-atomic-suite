      module scalve_i
      implicit none (type, external)
      interface
      subroutine scalve (vector, factor, ndim)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ndim
      real(real64), intent(inout) :: vector(ndim)
      real(real64), intent(in) :: factor
      end subroutine
      end interface
      end module
