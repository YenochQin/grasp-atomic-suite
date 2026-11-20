      module vecsum_i
      implicit none (type, external)
      interface
      subroutine vecsum (c, a, b, faca, facb, ndim)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ndim
      real(real64), intent(in)  :: a(ndim)
      real(real64), intent(in)  :: b(ndim)
      real(real64), intent(out) :: c(ndim)
      real(real64), intent(in) :: faca
      real(real64), intent(in) :: facb
      end subroutine vecsum
      end interface
      end module vecsum_i
