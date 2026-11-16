      module sbdat1_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,   only: max_orbitals
      implicit none (type, external)
!...Created by Pacific-Sierra Research 77to90  4.3E  06:25:32  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, parameter :: nlmax = 20
      integer :: nshlp(nlmax, nlmax)
      integer :: nshlpp(nlmax, max_orbitals)
      end module sbdat1_c
