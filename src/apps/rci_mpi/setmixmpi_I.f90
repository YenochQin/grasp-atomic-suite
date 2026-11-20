      module setmixmpi_i
      implicit none (type, external)
      interface
      subroutine setmixmpi (name, idblk)
      implicit none (type, external)
      character (len = *), intent(in) :: name
      character (len = 8) :: idblk(*)
      end subroutine setmixmpi
      end interface
      end module setmixmpi_i
