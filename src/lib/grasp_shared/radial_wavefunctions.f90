      module radial_wavefunctions
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only:  max_orbitals
      implicit none (type, external)
      private :: real64
      integer :: orbital_last_grid_index(max_orbitals)
      real(real64) :: orbital_origin_parameters(max_orbitals)
      integer :: mfff(max_orbitals)
      real(real64) :: pzff(max_orbitals)
      integer :: mfii(max_orbitals)
      real(real64) :: pzii(max_orbitals)
      ! Radial components are indexed by grid point, then by orbital.
      real(real64), pointer :: radial_large_components(:, :),radial_small_components(:, :)
      real(real64), pointer :: pfff(:, :),qfff(:, :),pfii(:, :),qfii(:, :)
      end module radial_wavefunctions
