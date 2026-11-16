      module sixj_i
      implicit none (type, external)
      interface
      subroutine sixj (two_j1, two_j2, two_j3, two_j4, two_j5, two_j6, check_triangle_conditions, six_j_value)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: two_j1
      integer, intent(in) :: two_j2
      integer, intent(in) :: two_j3
      integer, intent(in) :: two_j4
      integer, intent(in) :: two_j5
      integer, intent(in) :: two_j6
      integer, intent(in) :: check_triangle_conditions
      real(real64), intent(out) :: six_j_value
!VAST.../CONSTS/ ZERO(IN), ONE(IN)
!VAST...Calls: IXJTIK, SIXJ5, SIXJ1, SIXJ35, SIXJ2, GRACAH1
!VAST...Calls: SIXJ3, SIXJ4
      end subroutine
      end interface
      end module
