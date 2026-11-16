      module lagr_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:37:37  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: nec
      real(real64), pointer :: ecv(:)
      integer, pointer :: iecc(:)
      end module lagr_c
