      module invt_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,    only: max_orbitals
      implicit none (type, external)
      logical :: noinvt(max_orbitals)
      end module invt_c
