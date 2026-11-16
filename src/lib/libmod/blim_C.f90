      module blim_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      integer :: iprerun, ncsfpre
      real(real64) :: coeffcut1, coeffcut2
      end module blim_c
