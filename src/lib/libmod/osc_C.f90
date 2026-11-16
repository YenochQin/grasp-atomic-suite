      module osc_c
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  07:26:50   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: nsdim
      real(double), dimension(:), pointer :: hb1, hb2, hc1, hc2, hm1, hm2
      integer, dimension(:), pointer :: jja, jjb
      integer :: lk, kk
      integer :: ntdim
      real(double), dimension(:), pointer :: xsldr, totc, totb
      integer, dimension(:), pointer :: isldr, isldr1
      integer :: nint, ninteg
      integer, dimension(:), pointer :: nptr, lab
      integer :: nkp
      integer, dimension(:), pointer :: kp
      logical, dimension(10) :: ltc
      end module osc_c
