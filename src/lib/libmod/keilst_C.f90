      module keilst_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: ndkea, nkei
      integer, pointer :: indkei(:)
      real(real64), pointer :: valkei(:)
      logical :: frstki
      end module keilst_c
