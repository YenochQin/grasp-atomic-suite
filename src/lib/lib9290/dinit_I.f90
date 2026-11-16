      module dinit_i
      implicit none (type, external)
      interface
      subroutine dinit (n, a, x, incx)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: n
      real(real64), intent(in) :: a
      real(real64), intent(out) :: x(*)
      integer, intent(in) :: incx
      end subroutine
      end interface
      end module
