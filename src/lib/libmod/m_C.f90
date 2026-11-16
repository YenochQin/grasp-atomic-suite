!
!***********************************************************************
!                                                                      *
      module m_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,   only:  max_orbitals
      implicit none (type, external)
      integer :: nq1(max_orbitals), nq2(max_orbitals)
      integer :: jjc1(max_orbitals), jjc2(max_orbitals)
      integer :: jjq1(3, max_orbitals), jjq2(3, max_orbitals)
      integer :: jlist(max_orbitals), klist(max_orbitals)
      integer :: npeel, ncore
      end module m_c
