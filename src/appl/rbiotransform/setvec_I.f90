      module setvec_i
      implicit none (type, external)
      interface
      subroutine setvec (vector, value, ndim)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ndim
      real(real64), intent(out) :: vector(ndim)
      real(real64), intent(in) :: value
      end subroutine
      end interface
      end module
