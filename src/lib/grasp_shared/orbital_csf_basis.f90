! Stores orbital quantum numbers, energies, and occupations for the CSF basis.
      module orbital_csf_basis
      use, intrinsic :: iso_fortran_env, only: real64
      use, intrinsic :: iso_fortran_env, only: int8
      use grasp_parameters,   only:  max_orbitals
      implicit none (type, external)
      private :: real64, int8

      character(len=2) :: orbital_symmetry_labels(max_orbitals)
      character(len=2) :: nhr(max_orbitals)
      real(real64) :: orbital_energy_parameters(max_orbitals), orbital_origin_exponents(max_orbitals), ped(max_orbitals)
      integer :: csf_count, orbital_count, ncfr, nwr
      ! Occupation storage is indexed by orbital, then by CSF.
      integer(int8), pointer :: csf_subshell_occupations(:, :) => null()
      integer, pointer :: iqar(:, :) => null()
      integer :: orbital_principal_n(max_orbitals), orbital_kappa(max_orbitals), npr(max_orbitals), nakr(max_orbitals)
      integer :: orbital_l(max_orbitals), orbital_two_j(max_orbitals)
      end module orbital_csf_basis
