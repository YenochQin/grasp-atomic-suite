      module screenpar_i
      implicit none (type, external)
      interface
      subroutine screenpar (ncore)
      use parameter_def, only: max_orbitals
      implicit none (type, external)
      integer, intent(in) :: ncore
      end subroutine
      end interface
      end module
