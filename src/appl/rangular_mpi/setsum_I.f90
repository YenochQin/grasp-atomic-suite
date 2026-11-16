      module setsum_i
      implicit none (type, external)
      interface
      subroutine setsum (fullname)
      implicit none (type, external)
      character (len = *), intent(in) :: fullname
      end subroutine
      end interface
      end module
