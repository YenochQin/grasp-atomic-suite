      module transition_matrix_workspace
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
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
      end module transition_matrix_workspace
