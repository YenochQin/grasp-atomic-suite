      module osc_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  07:26:50   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: nsdim
      real(real64), pointer :: hb1(:), hb2(:), hc1(:), hc2(:), hm1(:), hm2(:)
      integer, pointer :: jja(:), jjb(:)
      integer :: lk, kk
      integer :: ntdim
      real(real64), pointer :: xsldr(:), totc(:), totb(:)
      integer, pointer :: isldr(:), isldr1(:)
      integer :: nint, ninteg
      integer, pointer :: nptr(:), lab(:)
      integer :: nkp
      integer, pointer :: kp(:)
      logical :: ltc(10)
      end module osc_c
