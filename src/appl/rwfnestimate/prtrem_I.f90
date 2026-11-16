      module prtrem_i
      implicit none (type, external)
      interface
      subroutine prtrem (all)
      use parameter_def, only: max_orbitals
      implicit none (type, external)
      logical, intent(out) :: all
      end subroutine
      end interface
      end module
