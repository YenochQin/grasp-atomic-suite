      module skint_i
      implicit none (type, external)
      interface
      real(real64) function skint (rac, ia, ic, rbd, ib, id, k, iw)
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_radial_points
      implicit none (type, external)
      real(real64), intent(in) :: rac(max_radial_points)
      integer, intent(in) :: ia
      integer, intent(in) :: ic
      real(real64), intent(in) :: rbd(max_radial_points)
      integer, intent(in) :: ib
      integer, intent(in) :: id
      integer, intent(in) :: k
      integer, intent(in) :: iw
      end function skint
      end interface
      end module skint_i
