      module setcsh_i
      implicit none (type, external)
      interface
      subroutine setcsh (nfile, name, ncore)
      implicit none (type, external)
      integer, intent(in) :: nfile
      character (len = *), intent(in) :: name
      integer, intent(out) :: ncore
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: LODCSH
!...This routine performs I/O.
      end subroutine
      end interface
      end module
