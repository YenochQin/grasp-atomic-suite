      module vinlst_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      integer :: ndvin, nvinti
      integer, pointer :: indtei(:)
      real(real64), pointer :: valtei(:)
      logical :: frstvi
      end module vinlst_c
