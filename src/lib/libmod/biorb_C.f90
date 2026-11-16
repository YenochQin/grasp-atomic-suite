      module biorb_c
      use vast_kind_param, only:double
      use parameter_def,   only: nnnw
!...Created by Pacific-Sierra Research 77to90  4.3E  06:22:09  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17

      character, dimension(nnnw) :: nhff*2, nhii*2
      real(double), dimension(nnnw) :: eff, gamaff , eii, gamaii
      integer :: ncfff, nwff, ncfii, nwii
      integer, dimension(nnnw) :: npff, nakff, npii, nakii
      integer, dimension(nnnw) :: nklff, nkjff, nklii, nkjii
      integer, dimension(nnnw) :: nnff, nnii
      end module biorb_c
