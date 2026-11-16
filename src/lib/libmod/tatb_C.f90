      module tatb_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,   only:  radial_workspace_capacity
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:23:52  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: mtp
      real(real64) :: ta(radial_workspace_capacity), tb(radial_workspace_capacity)
      end module tatb_c
