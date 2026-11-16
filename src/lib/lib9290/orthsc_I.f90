      module orthsc_i
      implicit none (type, external)
      interface
      subroutine orthsc
      implicit none (type, external)
!VAST.../DEBUGR/ LDBPR(IN)
!VAST.../DEF4/ ACCY(IN)
!VAST.../ORB2/ NW(IN)
!VAST.../ORB4/ NP(IN), NAK(IN)
!VAST.../ORB10/ NH(IN)
!VAST.../WAVE/ PZ(INOUT), MF(INOUT)
!VAST...Calls: RINT
!...This routine performs I/O.
      end subroutine
      end interface
      end module
