      module lodrwf_i
      implicit none (type, external)
      interface
      subroutine lodrwf (iierr)
      implicit none (type, external)
      integer, intent(out) :: iierr
!VAST.../DEBUGR/ LDBPR(IN)
!VAST.../DEF1/ Z(IN)
!VAST.../DEF2/ C(IN)
!VAST.../GRID/ N(IN)
!VAST.../NPAR/ NPARM(IN)
!VAST.../ORB1/ E(INOUT), GAMA(OUT)
!VAST.../ORB2/ NW(IN)
!VAST.../ORB4/ NP(IN), NAK(IN)
!VAST.../ORB10/ NH(IN)
!VAST.../WAVE/ PZ(OUT), PNTRPF(INOUT), PNTRQF(INOUT)
!VAST...Calls: ALLOC, INTRPQ, DALLOC, ORTHSC
!...This routine performs I/O.
      end subroutine
      end interface
      end module
