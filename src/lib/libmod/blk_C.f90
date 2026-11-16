!
!***********************************************************************
!                                                                      *
      module blk_c
!                                                                      *
!***********************************************************************
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  14:35:02   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, dimension(20) :: ncfblk, nevinblk, ncfinblk, two_j
      integer :: nblock
      real(double), dimension(:), pointer :: idxblk
      integer :: nelectot, ncftot, nwtot, nvectot, nvecsiztot, nblock1
      integer :: nblocki
      integer :: nblockf
      integer, dimension(10) :: ncfi
      integer, dimension(10) :: ncff
      end module blk_c
