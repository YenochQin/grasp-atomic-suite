      module frmhyd_i
      implicit none (type, external)
      interface
      subroutine frmhyd (matched_orbital_indices, matched_subshell_count, modify, nradial)
      use grasp_parameters, only: max_orbitals
      implicit none (type, external)
      integer, intent(in) :: matched_orbital_indices(max_orbitals)
      integer, intent(in) :: matched_subshell_count, nradial
      logical, intent(in) :: modify
      end subroutine frmhyd
      end interface
      end module frmhyd_i
