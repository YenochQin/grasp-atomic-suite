      module vinti_integral_table
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      integer :: ndvin, nvinti
      integer, pointer :: indtei(:)
      real(real64), pointer :: valtei(:)
      logical :: frstvi
      end module vinti_integral_table
