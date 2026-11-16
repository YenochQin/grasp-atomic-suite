      module setrwfmpi_i
      implicit none (type, external)
      interface
      subroutine setrwfmpi (name)
      implicit none (type, external)
      character :: name*(*)
      end subroutine
      end interface
      end module
