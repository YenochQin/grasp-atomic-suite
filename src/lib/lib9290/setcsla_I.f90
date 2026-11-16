      module setcsla_i
      implicit none (type, external)
      interface
      subroutine setcsla (name, ncore)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      integer, intent(out) :: ncore
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: OPENFL, LODCSL
!...This routine performs I/O.
      end subroutine
      end interface
      end module
