      module biorb_c
      use vast_kind_param, only:double
      use parameter_def,   only: max_orbitals
!...Created by Pacific-Sierra Research 77to90  4.3E  06:22:09  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17

      character, dimension(max_orbitals) :: nhff*2, nhii*2
      real(double), dimension(max_orbitals) :: eff, gamaff , eii, gamaii
      integer :: ncfff, nwff, ncfii, nwii
      integer, dimension(max_orbitals) :: npff, nakff, npii, nakii
      integer, dimension(max_orbitals) :: nklff, nkjff, nklii, nkjii
      integer, dimension(max_orbitals) :: nnff, nnii
      end module biorb_c
