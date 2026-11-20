!***********************************************************************
!                                                                      *
      subroutine lodmixmpi (block_labels)
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
      use atomic_calculation_state, only: electron_count, iccmin, ncmin
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use transition_basis_blocks,           only: idxblk
      use orbital_csf_basis,           only: csf_count, orbital_count, csf_subshell_occupations
      use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_integer4, mpi_rank
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
      integer :: nelec4, ncftot4, nw4, ncmin4, nvecsiz4, nblock4
!-----------------------------------------------

! lodstate generates
!    nevblk(), ncmaxblk()
!    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)
! Thus we let node-0 do it and then broadcast here

      call alloc (ncmaxblk, symmetry_block_count, "NCMAXBLK", "LODMIXmpi")
      call alloc (eigenvector_counts_by_block, symmetry_block_count, "NEVBLK", "LODMIXmpi")

!      print *, ' LODMIX: change lodstate arg-list - see rscf2/getold.f90 '
!      stop
      if (mpi_rank == 0) then
         call lodstate (block_labels)
      end if

      call mpi_bcast (eigenvector_counts_by_block, symmetry_block_count, mpi_integer, 0,               &
                                       mpi_comm_world, mpi_error_code)
      call mpi_bcast (ncmaxblk, symmetry_block_count, mpi_integer, 0,             &
                                       mpi_comm_world, mpi_error_code)
      call mpi_bcast (ncmin, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)

      if (ncmin /= 0) then
         if (mpi_rank /= 0) then   ! lodstate allocated it from node-0
            call alloc (iccmin, ncmin, "ICCMIN", "LODMIXmpi")
         end if

         call mpi_bcast (iccmin, ncmin, mpi_integer4, 0,            &
                                       mpi_comm_world, mpi_error_code)
      end if

!
!   Determine other auxiliary arrays, parameters
!
!    idxblk() is the block number of an eigenstate
!    nvecsiz  is the total size of the eigenvector array
!
      call alloc (idxblk, ncmin, "IDXBLK", "LODMIXmpi")
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
      if (noffset /= ncmin) stop "lodmixmpi: ncmin trouble"
!
!   Write header data to the  .mix  file
!
      if (mpi_rank == 0) then
         nelec4 = electron_count
         ncftot4 = total_csf_count
         nw4 = orbital_count
         ncmin4 = ncmin
         nvecsiz4 = vector_size
         nblock4 = symmetry_block_count
         write (25) nelec4, ncftot4, nw4, ncmin4, nvecsiz4, nblock4
      end if
!
      return
      end subroutine lodmixmpi
