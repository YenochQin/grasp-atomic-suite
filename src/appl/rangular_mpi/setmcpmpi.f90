!***********************************************************************
!                                                                      *
!cjb  myid, nprocs = NOT args
      subroutine setmcpmpi(ncore, idblk, filehead)
!                                                                      *
! A wrapper for setmcp/getinf. setmcp/getinf are then shared by serial *
! and MPI programs.                                                    *
!                                                                      *
!   Written by Xinghong He                Last revision: 30 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:01:42   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!cjb mpi_C
!     USE mpi_C,    ONLY: MYID, NPROCS, ierr
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use iccu_c
      use mcp_c
      use orb_c
      use def_c
      use hblock_c
!cjb mpi_C
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setmcp_i
!     USE setmcpmpi_I
      use getinf_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!      INTEGER , INTENT(IN) :: MYID
!      INTEGER , INTENT(IN) :: NPROCS
       integer , intent(in) :: ncore
      character , intent(in) :: filehead*(*)
       character , intent(in) :: idblk(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!cjb mpi_C
!cjb  INTEGER :: K, ierr
      integer :: k
!-----------------------------------------------

      call setmcp (ncore, idblk, filehead)

! DIAG, ICCUT, LFORDR are set in GETINF

      if (myid .eq. 0) then
         call getinf
      endif

      call mpi_bcast (iccut(1),100,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast (diag,  1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast (lfordr,1,mpi_logical,0,mpi_comm_world,ierr)

      do k = 30, 32+kmax
         write (k) nelec,ncf,nw
         write (k) diag,iccut(1),lfordr
      enddo

      do k = 30, 32+kmax
         close (k)
      enddo

      return
      end subroutine setmcpmpi
