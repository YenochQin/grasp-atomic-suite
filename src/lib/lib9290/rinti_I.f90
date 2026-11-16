      module rinti_i
      implicit none (type, external)
      interface
      real(real64) function rinti (j, k, mode)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: mode
!VAST.../DEBUGR/ LDBPR(IN)
!VAST.../DEF2/ C(IN)
!VAST.../GRID/ RP(IN), RPOR(IN), H(IN)
!VAST.../NPOT/ ZZ(IN)
!VAST.../ORB4/ NP(IN), NAK(IN)
!VAST.../ORB10/ NH(IN)
!VAST.../TATB/ TA(INOUT), TB(IN), MTP(OUT)
!VAST.../WAVE/ MF(IN)
!VAST...Calls: DPBDT, QF, PF, QUAD
!...This routine performs I/O.
      end function
      end interface
      end module
