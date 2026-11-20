!***********************************************************************
!                                                                      *
      subroutine lodmix (block_labels)
      use calculation_file_units, only: mixing_output_unit
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
      use array_memory_management, only: alloc
      use atomic_calculation_state, only: electron_count, ncmin
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use transition_basis_blocks, only: idxblk
      use orbital_csf_basis, only: csf_count, orbital_count, csf_subshell_occupations
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use block_level_selection, only: lodstate
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=8), intent(in):: block_labels(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: total_csf_count, noffset, vector_size, symmetry_block_index, level_index
!-----------------------------------------------

! lodstate generates
!    nevblk(), ncmaxblk()
!    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)
! Thus we let node-0 do it and then broadcast here

      call alloc (ncmaxblk, symmetry_block_count, "NCMAXBLK", "LODMIX")
      call alloc (eigenvector_counts_by_block, symmetry_block_count, "NEVBLK", "LODMIX")

!      print *, ' LODMIX: change lodstate arg-list - see rscf2/getold.f90 '
!      stop
      call lodstate (block_labels)
!     CALL lodstate (nblock, ncfblk(1), idblk, nevblk, ncmaxblk)
!
!   Determine other auxiliary arrays, parameters
!
!    idxblk() is the block number of an eigenstate
!    nvecsiz  is the total size of the eigenvector array
!
      call alloc (idxblk, ncmin, "IDXBLK", "LODMIX")
      total_csf_count = 0
      noffset = 0
      vector_size = 0
      do symmetry_block_index = 1, symmetry_block_count
         do level_index = 1, eigenvector_counts_by_block(symmetry_block_index)
            idxblk(level_index + noffset) = symmetry_block_index
         end do
         total_csf_count = total_csf_count + csf_counts_by_block(symmetry_block_index)
         noffset = noffset + eigenvector_counts_by_block(symmetry_block_index)
         vector_size = vector_size + csf_counts_by_block(symmetry_block_index) * eigenvector_counts_by_block(symmetry_block_index)
      end do
      if (noffset /= ncmin) stop "lodmix: ncmin trouble"
!
!   Write header data to the  .mix  file
!
      write (mixing_output_unit) electron_count, total_csf_count, orbital_count, ncmin, vector_size, symmetry_block_count
!
      return
      end subroutine lodmix
