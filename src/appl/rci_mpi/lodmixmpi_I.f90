      module lodmixmpi_i
      implicit none (type, external)
      interface
      subroutine lodmixmpi (idblk)
      implicit none (type, external)
      character(len=8), intent(in):: idblk(*)
      end subroutine lodmixmpi
      end interface
      end module lodmixmpi_i
