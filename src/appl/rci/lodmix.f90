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
      use blk_c, only: idxblk
      use orb_c, only: ncf, nw, iqa
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
!-----------------------------------------------

! lodstate generates
!    nevblk(), ncmaxblk()
!    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)
! Thus we let node-0 do it and then broadcast here

      call alloc (ncmaxblk, nblock, 'NCMAXBLK', 'LODMIX')
      call alloc (nevblk, nblock, 'NEVBLK', 'LODMIX')

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
      call alloc (idxblk, ncmin, 'IDXBLK', 'LODMIX')
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
      if (noffset .ne. ncmin) stop 'lodmix: ncmin trouble'
!
!   Write header data to the  .mix  file
!
      write (25) nelec, ncftot, nw, ncmin, nvecsiz, nblock
!
      return
      end subroutine lodmix
