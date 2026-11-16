      module tatb_c
      use vast_kind_param, only:  double
      use parameter_def,   only:  radial_workspace_capacity
!...Created by Pacific-Sierra Research 77to90  4.3E  06:23:52  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: mtp
      real(double), dimension(radial_workspace_capacity) :: ta, tb
      end module tatb_c
