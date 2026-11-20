      module setsum_i
      implicit none (type, external)
      interface
      subroutine setsum (name, nci)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      integer, intent(in) :: nci
      end subroutine setsum
      end interface
      end module setsum_i
