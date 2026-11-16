      module setsum_i
      implicit none (type, external)
      interface
      subroutine setsum (filnam)
      implicit none (type, external)
      character (len = *), intent(in) :: filnam
      end subroutine setsum
      end interface
      end module setsum_i
