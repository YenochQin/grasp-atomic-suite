      module orbout_i
      implicit none (type, external)
      interface
      subroutine orbout (rwffile2)
      implicit none (type, external)
      character (len = *), intent(in) :: rwffile2
      end subroutine
      end interface
      end module
