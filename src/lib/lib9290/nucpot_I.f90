      module nucpot_i
      implicit none (type, external)
      interface
      subroutine nucpot
      implicit none (type, external)
!VAST.../DEBUGR/ LDBPR(IN)
!VAST.../DEF1/ Z(IN)
!VAST.../DEF9/ PI(IN)
!VAST.../GRID/ R(IN), N(IN)
!VAST.../NPAR/ PARM(IN), NPARM(IN)
!VAST.../NPOT/ ZZ(INOUT), NNUC(OUT)
!VAST...Calls: ES
!...This routine performs I/O.
      end subroutine
      end interface
      end module
