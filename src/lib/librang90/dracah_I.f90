      module dracah_i
      implicit none (type, external)
      interface
      subroutine dracah (i, j, k, l, m, n, rac)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, parameter :: mfact = 500
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: l
      integer, intent(in) :: m
      integer, intent(in) :: n
      real(real64), intent(out) :: rac
!VAST.../FACTS/ GAM(IN)
      end subroutine
      end interface
      end module
