      module getrsl_i
      implicit none (type, external)
      interface
      subroutine getrsl (indx, nsubs)
      implicit none (type, external)
      integer, intent(inout) :: indx(*)
      integer, intent(out) :: nsubs
!VAST.../ORB2/ NW(IN)
!VAST.../ORB4/ NP(IN)
!VAST.../ORB10/ NH(IN)
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: LDIGIT, CONVRT
!...This routine performs I/O.
      end subroutine
      end interface
      end module
