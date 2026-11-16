!
!***********************************************************************
!                                                                      *
      module blk_c
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      integer :: ncfblk(20), nevinblk(20), ncfinblk(20), two_j(20)
      integer :: nblock
      real(real64), pointer :: idxblk(:)
      integer :: nelectot, ncftot, nwtot, nvectot, nvecsiztot, nblock1
      integer :: nblocki
      integer :: nblockf
      integer :: ncfi(10)
      integer :: ncff(10)
      end module blk_c
