      module invmat_i
      implicit none (type, external)
      interface
      subroutine invmat (a, b, matdim, ndim)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: matdim
      real(real64), intent(inout) :: a(matdim, matdim)
      real(real64) :: b(matdim, matdim)
      integer, intent(in) :: ndim
      end subroutine
      end interface
      end module
