      module frmhyd_i
      implicit none (type, external)
      interface
      subroutine frmhyd (index, nsubs, modify, nradial)
      use parameter_def, only: max_orbitals
      implicit none (type, external)
      integer, intent(in) :: index(max_orbitals)
      integer, intent(in) :: nsubs, nradial
      logical, intent(in) :: modify
      end subroutine frmhyd
      end interface
      end module frmhyd_i
