      module spodmv_i
      implicit none (type, external)
      interface
      subroutine spodmv (n, m, b, c)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(real64), intent(in) :: b(n, m)
      real(real64), intent(inout) :: c(n, m)
      end subroutine
      end interface
      end module
