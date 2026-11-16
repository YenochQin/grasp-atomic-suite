      module bess_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,   only:  max_radial_points
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(real64) :: wij(2)
      real(real64) :: bessj(2, 2, max_radial_points), bessn(2, 2, max_radial_points)
      real(real64) :: bj(max_radial_points, 3)
      real(real64) :: tc(max_radial_points), td(max_radial_points)
      end module bess_c
