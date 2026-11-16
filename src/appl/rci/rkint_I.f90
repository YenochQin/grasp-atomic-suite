      module rkint_i
      implicit none (type, external)
      interface
      real(real64) function rkint (rac, ia, ic, rbd, ib, id, k, iw)
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_radial_points
      implicit none (type, external)
      real(real64), intent(in) :: rac(max_radial_points)
      integer, intent(in) :: ia
      integer, intent(in) :: ic
      real(real64), intent(in) :: rbd(max_radial_points)
      integer, intent(in) :: ib
      integer, intent(in) :: id
      integer, intent(in) :: k
      integer, intent(in) :: iw
      end function rkint
      end interface
      end module rkint_i
