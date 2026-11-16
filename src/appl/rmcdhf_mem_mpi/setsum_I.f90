      module setsum_i
      implicit none (type, external)
      interface
      subroutine setsum (filnam)
      implicit none (type, external)
      character (len = *), intent(in) :: filnam
      end subroutine
      end interface
      end module
