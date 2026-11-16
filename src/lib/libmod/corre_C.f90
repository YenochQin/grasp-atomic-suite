      module corre_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,    only: max_orbitals
      implicit none (type, external)
      logical :: lcorre(max_orbitals)
      end module corre_c
