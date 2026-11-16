      module setpot_i
      implicit none (type, external)
      interface
      subroutine setpot (j, jp)
      implicit none (type, external)
      integer, intent(in) :: j
      integer, intent(out) :: jp
!VAST.../IOUNIT/ ISTDE(IN)
!VAST.../DEF2/ C(IN)
!VAST.../DEF4/ ACCY(IN)
!VAST.../GRID/ RP(IN), RPOR(IN), H(IN), N(IN)
!VAST.../INT3/ TF(OUT), TG(INOUT)
!VAST.../ORB1/ E(IN)
!VAST.../POTE/ YP(IN)
!...This routine performs I/O.
      end subroutine
      end interface
      end module
