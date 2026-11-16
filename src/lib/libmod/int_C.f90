      module int_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,    only: max_radial_points
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:37:37  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: mtp0
      real(real64) :: p(max_radial_points), q(max_radial_points)
      real(real64) :: p0, q0
      real(real64) :: tf(max_radial_points), tg(max_radial_points), xu(max_radial_points), xv(max_radial_points)
      end module int_c
