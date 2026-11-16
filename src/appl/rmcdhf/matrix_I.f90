      module matrix_i
      implicit none (type, external)
      interface
      subroutine matrix (dvdfirst)
      implicit none (type, external)
      logical, intent(in) ::  dvdfirst
      end subroutine matrix
      end interface
      end module matrix_i
