      module orba_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,    only: max_orbitals
      implicit none (type, external)
      integer :: iorder(max_orbitals)
      end module orba_c
