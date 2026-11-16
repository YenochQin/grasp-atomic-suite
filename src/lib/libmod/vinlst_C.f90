      module vinlst_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:35:13  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: ndvin, nvinti
      integer, pointer :: indtei(:)
      real(real64), pointer :: valtei(:)
      logical :: frstvi
      end module vinlst_c
