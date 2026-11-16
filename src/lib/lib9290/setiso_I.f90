      module setiso_i
      implicit none (type, external)
      interface
      subroutine setiso (fname)
      implicit none (type, external)
      character (len = *), intent(in) :: fname
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: OPENFL, LODISO
!...This routine performs I/O.
      end subroutine
      end interface
      end module
