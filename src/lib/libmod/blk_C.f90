!
!***********************************************************************
!                                                                      *
      module blk_c
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  14:35:02   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: ncfblk(20), nevinblk(20), ncfinblk(20), two_j(20)
      integer :: nblock
      real(real64), pointer :: idxblk(:)
      integer :: nelectot, ncftot, nwtot, nvectot, nvecsiztot, nblock1
      integer :: nblocki
      integer :: nblockf
      integer :: ncfi(10)
      integer :: ncff(10)
      end module blk_c
