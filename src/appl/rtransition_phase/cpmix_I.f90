      module cpmix_i
      implicit none (type, external)
      interface
      subroutine cpmix (name, inpci)
      implicit none (type, external)
      character (len = 24), intent(inout) :: name(2)
      integer, intent(in) :: inpci
      end subroutine
      end interface
      end module
