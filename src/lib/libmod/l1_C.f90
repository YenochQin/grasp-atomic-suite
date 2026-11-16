      module l1_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,   only:  max_orbitals
      implicit none (type, external)
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: jbq1(3, max_orbitals), jbq2(3, max_orbitals)
      integer :: jtq1(3), jtq2(3)
      end module l1_c
