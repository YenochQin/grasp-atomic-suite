      module bcore_c
!                     Gediminas Gaigalas  10/05/17
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,   only:  max_orbitals
      implicit none (type, external)
      integer :: icore(max_orbitals)
      end module bcore_c
