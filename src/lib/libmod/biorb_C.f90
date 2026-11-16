      module biorb_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,   only: max_orbitals
      implicit none (type, external)
      private :: real64

      character :: nhff(max_orbitals)*2, nhii(max_orbitals)*2
      real(real64) :: eff(max_orbitals), gamaff(max_orbitals) , eii(max_orbitals), gamaii(max_orbitals)
      integer :: ncfff, nwff, ncfii, nwii
      integer :: npff(max_orbitals), nakff(max_orbitals), npii(max_orbitals), nakii(max_orbitals)
      integer :: nklff(max_orbitals), nkjff(max_orbitals), nklii(max_orbitals), nkjii(max_orbitals)
      integer :: nnff(max_orbitals), nnii(max_orbitals)
      end module biorb_c
