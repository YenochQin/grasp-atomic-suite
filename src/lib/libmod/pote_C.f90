      module pote_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,    only: max_radial_points
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  11:02:52   1/ 2/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(real64) :: yp(max_radial_points), xp(max_radial_points), xq(max_radial_points)
      end module pote_c
