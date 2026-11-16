      module ltab_i
      implicit none (type, external)
      interface
      subroutine ltab (is, nqs, ks, irows)
      implicit none (type, external)
      integer, intent(in) :: is(4)
      integer, intent(inout) :: nqs(4)
      integer, intent(in) :: ks(4)
      integer, intent(inout) :: irows(4)
!VAST.../TERMS/ NROWS(IN)
!...This routine performs I/O.
      end subroutine
      end interface
      end module
