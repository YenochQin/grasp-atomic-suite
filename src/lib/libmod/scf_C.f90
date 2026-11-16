      module scf_c
      use vast_kind_param,  only: double
      use parameter_def,    only: max_orbitals
!...Created by Pacific-Sierra Research 77to90  4.3E  06:38:40  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(double), dimension(max_orbitals) :: ucf
      integer, dimension(max_orbitals) :: method
      real(double), dimension(max_orbitals) :: scnsty
      real(double) :: epsmin, epsmax, emin, emax, zinf
      integer :: ndcof, nxcof, nycof, nddim, nxdim, nydim
      real(double), dimension(:), pointer :: da, xa, ya
      integer, dimension(:), pointer :: nda, nxa, nya
      end module scf_c
