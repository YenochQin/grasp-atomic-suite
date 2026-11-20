      module frmtfp_i
      implicit none (type, external)
      interface
      subroutine frmtfp (matched_orbital_indices, matched_subshell_count)
      use grasp_parameters, only: max_orbitals
      implicit none (type, external)
      integer, intent(in) :: matched_orbital_indices(max_orbitals)
      integer, intent(in) :: matched_subshell_count
      end subroutine frmtfp
      end interface
      end module frmtfp_i
