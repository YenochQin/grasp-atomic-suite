      module radfile_i
      implicit none (type, external)
      interface
      subroutine radfile (name)
      implicit none (type, external)
      character (len = 24), intent(in) :: name(2)
      end subroutine radfile
      end interface
      end module radfile_i
