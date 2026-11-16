      module mohr_i
      implicit none (type, external)
      interface
      subroutine mohr (n, kappa, z, fzalfa)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, parameter :: numval = 12
      integer, intent(in) :: n
      integer, intent(in) :: kappa
      real(real64) :: z
      real(real64), intent(out) :: fzalfa
      end subroutine
      end interface
      end module
