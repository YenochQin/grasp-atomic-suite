      module frmrwf_i
      implicit none (type, external)
      interface
      subroutine frmrwf (matched_orbital_indices, matched_subshell_count, file_path)
      use grasp_parameters, only: radial_workspace_capacity, max_orbitals
      implicit none (type, external)
      integer, intent(in) :: matched_orbital_indices(max_orbitals)
      integer, intent(in) :: matched_subshell_count
      character (len = *), intent(in) :: file_path
      end subroutine frmrwf
      end interface
      end module frmrwf_i
