      module lodmixmpi_i
      implicit none (type, external)
      interface
      subroutine lodmixmpi (block_labels)
      implicit none (type, external)
      character(len=8), intent(in):: block_labels(*)
      end subroutine lodmixmpi
      end interface
      end module lodmixmpi_i
