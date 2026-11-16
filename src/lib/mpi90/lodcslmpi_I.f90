      module lodcslmpi_i
      implicit none (type, external)
      interface
      subroutine lodcslmpi (nfile, ncore, jblock)
      implicit none (type, external)
      integer :: nfile, ncore, jblock
      end subroutine
      end interface
      end module
