      module multbc_i
      implicit none (type, external)
      interface
      subroutine multbc (n, k, m, c, temp, b)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: n
      integer, intent(in) :: k
      integer, intent(in) :: m
      real(real64) :: c(k*m)
      real(real64) :: temp(m)
      real(real64) :: b(n*k)
!VAST...Calls: DGEMV, DCOPY
      end subroutine
      end interface
      end module
