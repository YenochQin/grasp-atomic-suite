      module coeils_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: ndcoea, ncoei
      integer, pointer :: indoei(:)
      real(real64), pointer :: valoei(:)
      logical :: frstco
      end module coeils_c
