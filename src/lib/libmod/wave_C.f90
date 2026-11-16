      module wave_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,   only:  max_orbitals
      implicit none (type, external)
      private :: real64
      integer :: orbital_last_grid_index(max_orbitals)
      real(real64) :: orbital_origin_parameters(max_orbitals)
      integer :: mfff(max_orbitals)
      real(real64) :: pzff(max_orbitals)
      integer :: mfii(max_orbitals)
      real(real64) :: pzii(max_orbitals)
      real(real64), pointer :: radial_large_components(:, :),radial_small_components(:, :)
      real(real64), pointer :: pfff(:, :),qfff(:, :),pfii(:, :),qfii(:, :)
      end module wave_c
