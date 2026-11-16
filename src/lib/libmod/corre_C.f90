      module corre_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,    only: max_orbitals
      implicit none (type, external)
!...Created by Pacific-Sierra Research 77to90  4.3E  15:25:01   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      logical :: lcorre(max_orbitals)
      end module corre_c
