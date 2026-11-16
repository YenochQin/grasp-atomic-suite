!***********************************************************************
!
      subroutine cslhmpi(name, ncore, nblkin, idblk)
!
!  A container which calls setcsll to open, read <name> file to get
!     nblock, ncfblk(), idblk(), ncftot.
!  It then calls lib92/lodcsh to get
!     ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
!  The file pointer points to the first CSL record after this routine.
!
!  Called by rcimpivu, mcpmpi
!
!  Xinghong He 98-06-23
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use memory_man
      use mpi_c
      use hblock_c
      use def_c, only: emn,ioncty,nelec,z
      use orb_c, ncftot=>ncf
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setcsll_i
      use lodcsh_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ncore
      integer, intent(in)  :: nblkin
      character  :: name*(*)
      character  :: idblk(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iqadum
!-----------------------------------------------
!

!     node-0 does exactly the same as the serial code does
      if (myid .eq. 0) then
        call alloc   (ncfblk, nblkin+1, 'ncfblk', 'cslhmpi')
        call setcsll (21, name, nblkin, nblock, ncfblk, ncftot, &
                       idblk)
        call ralloc  (ncfblk, nblock+1, 'ncfblk', 'cslhmpi')
         rewind (21)
         read   (21,*)
         !..Load header of <name> file
         call lodcsh (21, ncore)
      endif

! Broadcast results to other nodes. ncfblk should be allocated
! on these nodes (with myid .ne. 0)

      call mpi_bcast (nblock, 1, mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast (ncftot, 1, mpi_integer,0,mpi_comm_world,ierr)

      if (myid /= 0) then
         call alloc (ncfblk, 1+nblock, 'ncfblk', 'cslhmpi')
      endif

      call mpi_bcast (ncfblk, 1+nblock, mpi_integer,0, &
                                              mpi_comm_world,ierr)
      call mpi_bcast (idblk, 8*nblock, mpi_character,0, &
                                              mpi_comm_world,ierr)
      call mpi_bcast (ncore, 1, mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast (nelec, 1, mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast (nw,    1, mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast (np,   nw, mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast (nak,  nw, mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast (nkl,  nw, mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast (nkj,  nw, mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast (nh, 2*nw, mpi_character,0,mpi_comm_world,ierr)

      return
      end subroutine cslhmpi
