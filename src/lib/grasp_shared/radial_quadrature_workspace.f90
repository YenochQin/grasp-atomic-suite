      module radial_quadrature_workspace
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only:  radial_workspace_capacity
      implicit none (type, external)
      private :: real64
      integer :: mtp
      real(real64) :: ta(radial_workspace_capacity), tb(radial_workspace_capacity)
      end module radial_quadrature_workspace
