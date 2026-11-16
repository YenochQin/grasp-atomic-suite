      module vpilst_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:35:13  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: ndvpa, nvpi
      integer, pointer :: indvpi(:)
      real(real64), pointer :: valvpi(:)
      logical :: frstvp
      end module vpilst_c
