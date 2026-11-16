      module hydpar_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,   only: max_orbitals
      implicit none (type, external)
      private :: real64
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!...Created by Pacific-Sierra Research 77to90  4.3E  11:02:52   1/ 2/07
      real(real64) :: sigma(max_orbitals)
      end module hydpar_c
