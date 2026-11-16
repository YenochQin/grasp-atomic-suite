      module sbdat1_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,   only: max_orbitals
      implicit none (type, external)
      integer, parameter :: nlmax = 20
      integer :: nshlp(nlmax, nlmax)
      integer :: nshlpp(nlmax, max_orbitals)
      end module sbdat1_c
