      module fixd_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,    only: max_orbitals
      implicit none (type, external)
      integer :: fixed_orbital_count
      logical :: is_orbital_fixed(max_orbitals)
      end module fixd_c
