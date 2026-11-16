      module items_i
      implicit none (type, external)
      interface
      subroutine items (ncmin, ncf, record, ierr)
      implicit none (type, external)
      integer, intent(inout) :: ncmin
      integer, intent(in) :: ncf
      character (len = *), intent(in) :: record
      integer, intent(out) :: ierr
!VAST.../DEF7/ PCCMIN(INOUT)
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: ALLOC, CONVRT, RALLOC
!...This routine performs I/O.
      end subroutine
      end interface
      end module
