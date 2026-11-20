      module itjpo_i
      implicit none (type, external)
      interface
      integer function itjpo (csf_index)
      implicit none (type, external)
      integer, intent(in) :: csf_index
!VAST.../ORB2/ NCF(IN)
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: JCUPA
!...This routine performs I/O.
      end function itjpo
      end interface
      end module itjpo_i
