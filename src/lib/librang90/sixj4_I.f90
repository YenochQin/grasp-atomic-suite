      module sixj4_i
      implicit none (type, external)
      interface
      subroutine sixj4 (jc, je, jd, jb, jf, itik, si)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: jc
      integer, intent(in) :: je
      integer, intent(in) :: jd
      integer, intent(in) :: jb
      integer, intent(in) :: jf
      integer, intent(in) :: itik
      real(real64), intent(out) :: si
!VAST.../CONSTS/ ZERO(IN), HALF(IN), ONE(IN), TWO(IN), THREE(IN)
!VAST.../CONSTS/ EPS(IN)
!VAST...Calls: IXJTIK, GRACAH1, SIXJ2, SIXJ3
      end subroutine
      end interface
      end module
