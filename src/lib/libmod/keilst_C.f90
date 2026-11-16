      module keilst_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      integer :: ndkea, nkei
      integer, pointer :: indkei(:)
      real(real64), pointer :: valkei(:)
      logical :: frstki
      end module keilst_c
