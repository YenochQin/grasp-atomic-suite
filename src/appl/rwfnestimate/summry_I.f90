      module summry_i
      implicit none (type, external)
      interface
      subroutine summry (nunit)
      use grasp_parameters, only: radial_workspace_capacity, max_orbitals
      implicit none (type, external)
      integer, intent(in) :: nunit
      end subroutine summry
      end interface
      end module summry_i
