      module wave_c
      use vast_kind_param, only:  double
      use parameter_def,   only:  nnnw
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, dimension(nnnw) :: mf
      real(double), dimension(nnnw) :: pz
      integer, dimension(nnnw) :: mfff
      real(double), dimension(nnnw) :: pzff
      integer, dimension(nnnw) :: mfii
      real(double), dimension(nnnw) :: pzii
      real(double), dimension(:,:), pointer :: pf,qf
      real(double), dimension(:,:), pointer :: pfff,qfff,pfii,qfii
      end module wave_c
