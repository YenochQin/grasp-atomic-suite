      module l1_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,   only:  max_orbitals
      implicit none (type, external)
      integer :: jbq1(3, max_orbitals), jbq2(3, max_orbitals)
      integer :: jtq1(3), jtq2(3)
      end module l1_c
