      module horb_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,   only: max_radial_points
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(real64) :: ph(max_radial_points), qh(max_radial_points)
      end module horb_c
