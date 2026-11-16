      module buffer_c
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  13:03:28   1/25/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: nbdim, nvcoef
      integer, dimension(:,:), pointer :: label
      real(double), dimension(:), pointer :: coeff
!     Active only in the audited rangular OpenMP object build.
!$OMP THREADPRIVATE(NBDIM, NVCOEF, LABEL, COEFF)
      end module buffer_c
