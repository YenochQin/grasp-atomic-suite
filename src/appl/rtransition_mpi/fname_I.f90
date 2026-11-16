      module fname_i
      implicit none (type, external)
      interface
      subroutine fname (name)
      implicit none (type, external)
      character (len = 24), intent(in) :: name(2)
!...This routine performs I/O.
      end subroutine fname
      end interface
      end module fname_i
