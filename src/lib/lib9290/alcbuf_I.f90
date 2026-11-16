      module alcbuf_i
      implicit none (type, external)
      interface
      subroutine alcbuf (mode)
      implicit none (type, external)
      integer, intent(in) :: mode
!VAST...Calls: ALLOC, RALLOC, DALLOC
!...This routine performs I/O.
      end subroutine
      end interface
      end module
