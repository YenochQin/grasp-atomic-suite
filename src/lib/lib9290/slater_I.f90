      module slater_i
      implicit none (type, external)
      interface
      real(real64) function slater (ia, ib, ic, id, k)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ia
      integer, intent(in) :: ib
      integer, intent(in) :: ic
      integer, intent(in) :: id
      integer, intent(in) :: k
!VAST.../DEBUGR/ LDBPR(IN)
!VAST.../GRID/ RPOR(IN)
!VAST.../ORB4/ NP(IN)
!VAST.../ORB10/ NH(IN)
!VAST.../TATB/ TA(OUT), TB(IN), MTP(INOUT)
!VAST.../WAVE/ MF(IN)
!VAST...Calls: YZK, PF, QF, QUAD
!...This routine performs I/O.
      end function
      end interface
      end module
