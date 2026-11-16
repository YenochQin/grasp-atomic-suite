      module wrtmat_i
      implicit none (type, external)
      interface
      subroutine wrtmat (a, nrow, ncol, nmrow, nmcol)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: nmrow
      integer, intent(in) :: nmcol
      real(real64), intent(in) :: a(nmrow, nmcol)
      integer, intent(in) :: nrow
      integer, intent(in) :: ncol
!...This routine performs I/O.
      end subroutine
      end interface
      end module
