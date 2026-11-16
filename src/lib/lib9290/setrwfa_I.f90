      module setrwfa_i
      implicit none (type, external)
      interface
      subroutine setrwfa (name)
      implicit none (type, external)
      character (len = *), intent(in) :: name
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: OPENFL, LODRWF
!...This routine performs I/O.
      end subroutine
      end interface
      end module
