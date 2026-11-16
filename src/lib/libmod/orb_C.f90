!
!***********************************************************************
!                                                                      *
      module orb_c
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64, double => real64, byte => int8
      use, intrinsic :: iso_fortran_env, only: int8
      use parameter_def,   only:  max_orbitals
      implicit none (type, external)
      private :: real64, int8

      character(len=2) :: orbital_symmetry_labels(max_orbitals)
      character(len=2) :: nhr(max_orbitals)
      real(real64) :: orbital_energy_parameters(max_orbitals), orbital_origin_exponents(max_orbitals), ped(max_orbitals)
      integer :: csf_count, orbital_count, ncfr, nwr
      integer(int8), pointer :: csf_subshell_occupations(:, :)
      integer, pointer :: iqar(:, :)
      integer :: orbital_principal_n(max_orbitals), orbital_kappa(max_orbitals), npr(max_orbitals), nakr(max_orbitals)
      integer :: orbital_l(max_orbitals), orbital_two_j(max_orbitals)
      end module orb_c
