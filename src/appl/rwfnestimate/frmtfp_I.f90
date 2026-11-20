      module frmtfp_i
      implicit none (type, external)
      interface
      subroutine frmtfp (index, nsubs)
      use grasp_parameters, only: max_orbitals
      implicit none (type, external)
      integer, intent(in) :: index(max_orbitals)
      integer, intent(in) :: nsubs
      end subroutine frmtfp
      end interface
      end module frmtfp_i
