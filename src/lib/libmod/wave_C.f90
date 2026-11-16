      module wave_c
      use vast_kind_param, only:  double
      use parameter_def,   only:  max_orbitals
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, dimension(max_orbitals) :: orbital_last_grid_index
      real(double), dimension(max_orbitals) :: orbital_origin_parameters
      integer, dimension(max_orbitals) :: mfff
      real(double), dimension(max_orbitals) :: pzff
      integer, dimension(max_orbitals) :: mfii
      real(double), dimension(max_orbitals) :: pzii
      real(double), dimension(:,:), pointer :: radial_large_components,radial_small_components
      real(double), dimension(:,:), pointer :: pfff,qfff,pfii,qfii
      end module wave_c
