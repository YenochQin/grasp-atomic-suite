      module setsum_i
      implicit none (type, external)
      interface
      subroutine setsum (fullname)
      implicit none (type, external)
      character (len = *), intent(in) :: fullname
      end subroutine setsum
      end interface
      end module setsum_i
