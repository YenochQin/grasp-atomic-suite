      module itjpo_i
      implicit none (type, external)
      interface
      integer function itjpo (icsf)
      implicit none (type, external)
      integer, intent(in) :: icsf
!VAST.../ORB2/ NCF(IN)
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: JCUPA
!...This routine performs I/O.
      end function itjpo
      end interface
      end module itjpo_i
