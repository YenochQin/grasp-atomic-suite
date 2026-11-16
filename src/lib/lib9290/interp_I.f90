      module interp_i
      implicit none (type, external)
      interface
      subroutine interp (xarr, yarr, narr, xval, yval, accy)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, parameter :: mxord = 11
      integer, intent(in) :: narr
      real(real64), intent(in) :: xarr(narr)
      real(real64), intent(in) :: yarr(narr)
      real(real64), intent(in) :: xval
      real(real64), intent(out) :: yval
      real(real64), intent(in) :: accy
!...This routine performs I/O.
      end subroutine
      end interface
      end module
