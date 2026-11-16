      module strsum_i
      implicit none (type, external)
      interface
      subroutine strsum (name, inpci,ilbl)
      implicit none (type, external)
      character (len = 24), intent(in) :: name(2)
      integer, intent(in) :: inpci
      integer :: ilbl
      end subroutine
      end interface
      end module
