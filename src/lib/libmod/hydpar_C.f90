      module hydpar_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,   only: max_orbitals
      implicit none (type, external)
      private :: real64
      real(real64) :: sigma(max_orbitals)
      end module hydpar_c
