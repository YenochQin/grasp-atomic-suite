      module hblock_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      implicit none (type, external)
      integer :: symmetry_block_count  !! NCFBLK as an integer needs to be renamed
!     REAL(DOUBLE) :: PNCFBLK
!     REAL(DOUBLE) :: PNEVBLK, PNCMAXBLK
      integer, pointer :: csf_counts_by_block(:) => null()  !! this is a problem
      integer, pointer :: eigenvector_counts_by_block(:), ncmaxblk(:)
      end module hblock_c
