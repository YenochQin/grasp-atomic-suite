!***********************************************************************
!                                                                      *
      subroutine lodmixmpi (idblk)
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
      use hblock_c
      use blk_c,           only: idxblk
      use orb_c,           only: csf_count, orbital_count, csf_subshell_occupations
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lodstate_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=8), intent(in):: idblk(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncftot, noffset, nvecsiz, jb, j
      integer :: nelec4, ncftot4, nw4, ncmin4, nvecsiz4, nblock4
!-----------------------------------------------

! lodstate generates
!    nevblk(), ncmaxblk()
!    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)
! Thus we let node-0 do it and then broadcast here

      call alloc (ncmaxblk, symmetry_block_count, 'NCMAXBLK', 'LODMIXmpi')
      call alloc (eigenvector_counts_by_block, symmetry_block_count, 'NEVBLK', 'LODMIXmpi')

!      print *, ' LODMIX: change lodstate arg-list - see rscf2/getold.f90 '
!      stop
      if (mpi_rank .eq. 0) then
         call lodstate (idblk)
      endif

      call mpi_bcast (eigenvector_counts_by_block, symmetry_block_count, mpi_integer, 0,               &
                                       mpi_comm_world, mpi_error_code)
      call mpi_bcast (ncmaxblk, symmetry_block_count, mpi_integer, 0,             &
                                       mpi_comm_world, mpi_error_code)
      call mpi_bcast (ncmin, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)

      if (ncmin .ne. 0) then
         if (mpi_rank .ne. 0) then   ! lodstate allocated it from node-0
            call alloc (iccmin, ncmin, 'ICCMIN', 'LODMIXmpi')
         endif

         call mpi_bcast (iccmin, ncmin, mpi_integer4, 0,            &
                                       mpi_comm_world, mpi_error_code)
      endif

!
!   Determine other auxiliary arrays, parameters
!
!    idxblk() is the block number of an eigenstate
!    nvecsiz  is the total size of the eigenvector array
!
      call alloc (idxblk, ncmin, 'IDXBLK', 'LODMIXmpi')
      ncftot = 0
      noffset = 0
      nvecsiz = 0
      do jb = 1, symmetry_block_count
         do j = 1, eigenvector_counts_by_block(jb)
            idxblk(j + noffset) = jb
         enddo
         ncftot = ncftot + csf_counts_by_block(jb)
         noffset = noffset + eigenvector_counts_by_block(jb)
         nvecsiz = nvecsiz + csf_counts_by_block(jb) * eigenvector_counts_by_block(jb)
      enddo
      if (noffset .ne. ncmin) stop 'lodmixmpi: ncmin trouble'
!
!   Write header data to the  .mix  file
!
      if (mpi_rank .eq. 0) then
         nelec4 = electron_count
         ncftot4 = ncftot
         nw4 = orbital_count
         ncmin4 = ncmin
         nvecsiz4 = nvecsiz
         nblock4 = symmetry_block_count
         write (25) nelec4, ncftot4, nw4, ncmin4, nvecsiz4, nblock4
      endif
!
      return
      end subroutine lodmixmpi
