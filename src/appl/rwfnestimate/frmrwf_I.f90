      module frmrwf_i
      implicit none (type, external)
      interface
      subroutine frmrwf (index, nsubs, filnam)
      use parameter_def, only: radial_workspace_capacity, max_orbitals
      implicit none (type, external)
      integer, intent(in) :: index(max_orbitals)
      integer, intent(in) :: nsubs
      character (len = *), intent(in) :: filnam
      end subroutine
      end interface
      end module
