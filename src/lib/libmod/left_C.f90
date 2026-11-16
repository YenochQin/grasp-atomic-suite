      module left_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,   only:  max_orbitals
      implicit none (type, external)
      logical :: set(max_orbitals)
      end module left_c
