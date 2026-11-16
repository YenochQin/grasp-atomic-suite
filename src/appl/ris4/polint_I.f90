      module polint_i
      implicit none (type, external)
      interface
      subroutine polint (xa, ya, dens)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: xa(3)
      real(real64), intent(in) :: ya(3)
      real(real64), intent(out) :: dens
      end subroutine polint
      end interface
      end module polint_i
