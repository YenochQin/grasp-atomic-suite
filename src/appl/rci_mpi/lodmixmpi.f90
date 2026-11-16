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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use memory_man
      use def_c
      use hblock_c
      use blk_c,           only: idxblk
      use orb_c,           only: ncf, nw, iqa
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lodstate_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=8),dimension(*), intent(in):: idblk
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

      call alloc (ncmaxblk, nblock, 'NCMAXBLK', 'LODMIXmpi')
      call alloc (nevblk, nblock, 'NEVBLK', 'LODMIXmpi')

!      print *, ' LODMIX: change lodstate arg-list - see rscf2/getold.f90 '
!      stop
      if (myid .eq. 0) then
         call lodstate (idblk)
      endif

      call mpi_bcast (nevblk, nblock, mpi_integer, 0,               &
                                       mpi_comm_world, ierr)
      call mpi_bcast (ncmaxblk, nblock, mpi_integer, 0,             &
                                       mpi_comm_world, ierr)
      call mpi_bcast (ncmin, 1, mpi_integer, 0, mpi_comm_world, ierr)

      if (ncmin .ne. 0) then
         if (myid .ne. 0) then   ! lodstate allocated it from node-0
            call alloc (iccmin, ncmin, 'ICCMIN', 'LODMIXmpi')
         endif

         call mpi_bcast (iccmin, ncmin, mpi_integer4, 0,            &
                                       mpi_comm_world, ierr)
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
      do jb = 1, nblock
         do j = 1, nevblk(jb)
            idxblk(j + noffset) = jb
         enddo
         ncftot = ncftot + ncfblk(jb)
         noffset = noffset + nevblk(jb)
         nvecsiz = nvecsiz + ncfblk(jb) * nevblk(jb)
      enddo
      if (noffset .ne. ncmin) stop 'lodmixmpi: ncmin trouble'
!
!   Write header data to the  .mix  file
!
      if (myid .eq. 0) then
         nelec4 = nelec
         ncftot4 = ncftot
         nw4 = nw
         ncmin4 = ncmin
         nvecsiz4 = nvecsiz
         nblock4 = nblock
         write (25) nelec4, ncftot4, nw4, ncmin4, nvecsiz4, nblock4
      endif
!
      return
      end subroutine lodmixmpi
