      module hblock_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      implicit none (type, external)
!...Created by Pacific-Sierra Research 77to90  4.3E  06:29:39  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: symmetry_block_count  !! NCFBLK as an integer needs to be renamed
!     REAL(DOUBLE) :: PNCFBLK
!     REAL(DOUBLE) :: PNEVBLK, PNCMAXBLK
      integer, pointer :: csf_counts_by_block(:)  !! this is a problem
      integer, pointer :: eigenvector_counts_by_block(:), ncmaxblk(:)
      end module hblock_c
