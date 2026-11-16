      module whfrom_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,   only: max_orbitals
      implicit none (type, external)
!...Created by Pacific-Sierra Research 77to90  4.3E  11:02:52   1/ 2/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      character :: source(max_orbitals)*256
      end module whfrom_c
