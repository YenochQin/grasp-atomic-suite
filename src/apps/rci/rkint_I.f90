      module rkint_i
      implicit none (type, external)
      interface
      real(real64) function rkint (recoupling_coefficient, orbital_a, orbital_c, rbd, orbital_b, orbital_d, multipole_rank, iw)
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_radial_points
      implicit none (type, external)
      real(real64), intent(in) :: recoupling_coefficient(max_radial_points)
      integer, intent(in) :: orbital_a
      integer, intent(in) :: orbital_c
      real(real64), intent(in) :: rbd(max_radial_points)
      integer, intent(in) :: orbital_b
      integer, intent(in) :: orbital_d
      integer, intent(in) :: multipole_rank
      integer, intent(in) :: iw
      end function rkint
      end interface
      end module rkint_i
