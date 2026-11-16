      module vpilst_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      integer :: ndvpa, nvpi
      integer, pointer :: indvpi(:)
      real(real64), pointer :: valvpi(:)
      logical :: frstvp
      end module vpilst_c
