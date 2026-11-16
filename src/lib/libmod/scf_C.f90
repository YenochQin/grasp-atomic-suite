      module scf_c
      use vast_kind_param,  only: double
      use parameter_def,    only: nnnw
!...Created by Pacific-Sierra Research 77to90  4.3E  06:38:40  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(double), dimension(nnnw) :: ucf
      integer, dimension(nnnw) :: method
      real(double), dimension(nnnw) :: scnsty
      real(double) :: epsmin, epsmax, emin, emax, zinf
      integer :: ndcof, nxcof, nycof, nddim, nxdim, nydim
      real(double), dimension(:), pointer :: da, xa, ya
      integer, dimension(:), pointer :: nda, nxa, nya
      end module scf_c
