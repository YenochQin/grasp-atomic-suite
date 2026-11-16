      module parsjl_i
      implicit none (type, external)
      interface
      subroutine parsjl (mode, ncore, record, loc, jx, njx, ierr)
      implicit none (type, external)
      integer, intent(in) :: mode
      integer, intent(in) :: ncore
      character (len = 256), intent(in) :: record
      integer, intent(in) :: loc
      integer, intent(out) :: jx(*)
      integer, intent(out) :: njx
      integer, intent(out) :: ierr
!VAST.../ORB2/ NW(IN)
!VAST...Calls: CONVRT
!...This routine performs I/O.
      end subroutine
      end interface
      end module
