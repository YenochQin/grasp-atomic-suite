! Supplies the logarithmic factorial table used by the angular-momentum algebra.
      module log_factorial_cache
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      integer, parameter :: mfact = 500
      real(real64) :: gam(mfact)
      end module log_factorial_cache
