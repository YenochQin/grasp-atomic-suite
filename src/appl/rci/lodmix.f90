!***********************************************************************
!                                                                      *
      subroutine lodmix (idblk)
!                                                                      *
!   Determines the eigenpairs required;  this information is written   *
!   to the head of the  .mix  file.                                    *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, lodstate                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 24 Nov 1992   *
!   Block version by Xinghong He          Last revision:  9 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use memory_man
      use def_c
      use hblock_c, only: csf_counts_by_block, double, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use blk_c, only: idxblk
      use orb_c, only: csf_count, orbital_count, csf_subshell_occupations
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lodstate_i, only: lodstate
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=8), intent(in):: idblk(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncftot, noffset, nvecsiz, jb, j
!-----------------------------------------------

! lodstate generates
!    nevblk(), ncmaxblk()
!    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)
! Thus we let node-0 do it and then broadcast here

      call alloc (ncmaxblk, symmetry_block_count, "NCMAXBLK", "LODMIX")
      call alloc (eigenvector_counts_by_block, symmetry_block_count, "NEVBLK", "LODMIX")

!      print *, ' LODMIX: change lodstate arg-list - see rscf2/getold.f90 '
!      stop
      call lodstate (idblk)
!     CALL lodstate (nblock, ncfblk(1), idblk, nevblk, ncmaxblk)
!
!   Determine other auxiliary arrays, parameters
!
!    idxblk() is the block number of an eigenstate
!    nvecsiz  is the total size of the eigenvector array
!
      call alloc (idxblk, ncmin, "IDXBLK", "LODMIX")
      ncftot = 0
      noffset = 0
      nvecsiz = 0
      do jb = 1, symmetry_block_count
         do j = 1, eigenvector_counts_by_block(jb)
            idxblk(j + noffset) = jb
         end do
         ncftot = ncftot + csf_counts_by_block(jb)
         noffset = noffset + eigenvector_counts_by_block(jb)
         nvecsiz = nvecsiz + csf_counts_by_block(jb) * eigenvector_counts_by_block(jb)
      end do
      if (noffset /= ncmin) stop "lodmix: ncmin trouble"
!
!   Write header data to the  .mix  file
!
      write (25) electron_count, ncftot, orbital_count, ncmin, nvecsiz, symmetry_block_count
!
      return
      end subroutine lodmix
