      module convrt2_i
      implicit none (type, external)
      interface
      subroutine convrt2 (intnum, cnum, lenth, from)
      implicit none (type, external)
      integer, intent(in) :: intnum
      character (len = *), intent(inout) :: cnum
      integer, intent(out) :: lenth
      character (len = *), intent(in) :: from
!VAST.../IOUNIT/ ISTDE(IN)
!...This routine performs I/O.
      end subroutine
      end interface
      end module
