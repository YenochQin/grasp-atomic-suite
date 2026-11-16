      module sixj2_i
      implicit none (type, external)
      interface
      subroutine sixj2 (j, k, l, m, n, itik, si)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: l
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer, intent(in) :: itik
      real(real64), intent(out) :: si
!VAST.../CONSTS/ ZERO(IN), HALF(IN), ONE(IN), TWO(IN), THREE(IN)
!VAST.../CONSTS/ FOUR(IN)
!VAST...Calls: IXJTIK
      end subroutine
      end interface
      end module
