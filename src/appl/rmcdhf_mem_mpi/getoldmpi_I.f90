      module getoldmpi_i
      implicit none (type, external)
      interface
      subroutine getoldmpi (idblk)
      implicit none (type, external)
      character (len = 8) :: idblk(*)
      end subroutine
      end interface
      end module
