      module bcore_c
!                     Gediminas Gaigalas  10/05/17
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,   only:  max_orbitals
      implicit none (type, external)
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
      integer :: icore(max_orbitals)
      end module bcore_c
