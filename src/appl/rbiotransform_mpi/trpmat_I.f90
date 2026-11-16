      module trpmat_i
      implicit none (type, external)
      interface
      subroutine trpmat (xin, nrow, ncol, xout)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: nrow
      integer, intent(in) :: ncol
      real(real64), intent(in) :: xin(nrow, ncol)
      real(real64), intent(out) :: xout(ncol, nrow)
      end subroutine
      end interface
      end module
