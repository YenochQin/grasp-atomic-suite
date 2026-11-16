      module grid_c
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  06:22:09  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      use parameter_def,   only: nnn1
      integer :: n
      real(double), dimension(nnn1) :: r, rp, rpor
      real(double) :: rnt, h, hp
      end module grid_c
