      module lodstate_i
      implicit none (type, external)
      interface
      subroutine lodstate(idblk)
      implicit none (type, external)
      character (len = 8) :: idblk(*)
!VAST.../DEF7/ PCCMIN(INOUT), NCMIN(OUT)
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: ITEMS, DALLOC, ICCMIN
!...This routine performs I/O.
      end subroutine
      end interface
      end module
