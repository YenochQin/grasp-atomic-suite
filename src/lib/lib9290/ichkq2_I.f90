      module ichkq2_i
      implicit none (type, external)
      interface
      integer function ichkq2 (ja, jb)
      implicit none (type, external)
      integer, intent(in) :: ja
      integer, intent(in) :: jb
!VAST.../ORB2/ NW(IN)
!VAST...Calls: IQ
      end function
      end interface
      end module
