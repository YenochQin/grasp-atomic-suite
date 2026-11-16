      module jcup_i
      implicit none (type, external)
      interface
      integer function jcup (loc, icsf)
      implicit none (type, external)
      integer, intent(in) :: loc
      integer, intent(in) :: icsf
!VAST.../ORB2/ NCF(IN), NW(IN)
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: JCUPA
!...This routine performs I/O.
      end function
      end interface
      end module
