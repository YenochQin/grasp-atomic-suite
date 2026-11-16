!
!***********************************************************************
!                                                                      *
      module orb_c
!                                                                      *
!***********************************************************************
      use vast_kind_param, only:  double, byte
      use parameter_def,   only:  max_orbitals
!...Created by Pacific-Sierra Research 77to90  4.3E  08:57:22  12/25/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17

      character(len=2), dimension(max_orbitals) :: orbital_symmetry_labels
      character(len=2), dimension(max_orbitals) :: nhr
      real(double), dimension(max_orbitals) :: orbital_energy_parameters, orbital_origin_exponents, ped
      integer :: csf_count, orbital_count, ncfr, nwr
      integer(byte), dimension(:,:), pointer :: csf_subshell_occupations
      integer, dimension(:,:), pointer :: iqar
      integer, dimension(max_orbitals) :: orbital_principal_n, orbital_kappa, npr, nakr
      integer, dimension(max_orbitals) :: orbital_l, orbital_two_j
      end module orb_c
