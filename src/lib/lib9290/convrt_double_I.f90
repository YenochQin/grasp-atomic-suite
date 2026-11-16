      module convrt_double_i
      implicit none (type, external)
      interface
      subroutine convrt_double (intnum, cnum, lenth)
      implicit none (type, external)
      integer, intent(in) :: intnum
      character (len = *), intent(inout) :: cnum
      integer, intent(out) :: lenth
!...This routine performs I/O.
      end subroutine
      end interface
      end module
