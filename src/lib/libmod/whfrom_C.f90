      module whfrom_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,   only: max_orbitals
      implicit none (type, external)
      character :: source(max_orbitals)*256
      end module whfrom_c
