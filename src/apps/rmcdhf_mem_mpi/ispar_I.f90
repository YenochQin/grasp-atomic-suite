      module ispar_i
      implicit none (type, external)
      interface
      integer function ispar (csf_index)
      implicit none (type, external)
      integer, intent(in) :: csf_index
!VAST.../ORB2/ NCF(IN)
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: JCUPA
!...This routine performs I/O.
      end function ispar
      end interface
      end module ispar_i
