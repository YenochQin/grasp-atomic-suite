      module hmat_c
      use vast_kind_param, only:  double, long
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer(long) :: nelmnt
      real(double), dimension(:), pointer :: emt
      integer, dimension(:), pointer :: iendc
      integer, dimension(:), pointer :: irow
      integer, dimension(6) :: ntpitmp
      integer(long) :: nelmnttmp
      integer :: ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp,  &
                 ncoretmp, nvpitmp, nkeitmp, nvintitmp, ncftmp
      real(double) :: cutofftmp
      end module hmat_c
