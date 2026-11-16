      module tatb_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,   only:  radial_workspace_capacity
      implicit none (type, external)
      private :: real64
      integer :: mtp
      real(real64) :: ta(radial_workspace_capacity), tb(radial_workspace_capacity)
      end module tatb_c
