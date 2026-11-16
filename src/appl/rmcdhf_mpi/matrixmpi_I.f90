      module matrixmpi_i
      implicit none (type, external)
      interface
      subroutine matrixmpi (dvdfirst)
      implicit none (type, external)
      logical, intent(in) ::  dvdfirst
      end subroutine matrixmpi
      end interface
      end module matrixmpi_i
