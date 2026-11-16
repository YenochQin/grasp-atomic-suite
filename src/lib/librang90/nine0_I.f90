      module nine0_i
      implicit none (type, external)
      interface
      subroutine nine0 (j1, j2, j3, l1, l2, l3, k1, k2, k3, aa)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: j1
      integer, intent(in) :: j2
      integer, intent(in) :: j3
      integer, intent(in) :: l1
      integer, intent(in) :: l2
      integer, intent(in) :: l3
      integer, intent(in) :: k1
      integer, intent(in) :: k2
      integer, intent(in) :: k3
      real(real64), intent(out) :: aa
!VAST.../CONSTS/ ZERO(IN), ONE(IN)
!VAST...Calls: SIXJ
      end subroutine
      end interface
      end module
