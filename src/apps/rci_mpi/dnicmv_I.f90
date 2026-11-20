      module dnicmv_i
      implicit none (type, external)
      interface
      subroutine dnicmv (n, m, b, c)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(real64), intent(in) :: b(n, m)
      real(real64), intent(inout) :: c(n, m)
      end subroutine dnicmv
      end interface
      end module dnicmv_i
