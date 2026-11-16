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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!cjb mpi_C
!     USE mpi_C,    ONLY: MYID, NPROCS, ierr
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use iccu_c, only: iccut, iccutblk
      use mcp_c, only: diagonal_only, double, max_coulomb_multipole_rank, use_first_order_restriction
      use orb_c
      use def_c
      use hblock_c, only: csf_counts_by_block, double, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
!cjb mpi_C
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setmcp_i, only: setmcp
!     USE setmcpmpi_I
      use getinf_i, only: getinf
      implicit none (type, external)
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

      if (mpi_rank == 0) then
         call getinf
      end if

      call mpi_bcast (iccut(1),100,mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (diagonal_only,  1,mpi_logical,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (use_first_order_restriction,1,mpi_logical,0,mpi_comm_world,mpi_error_code)

      do k = 30, 32+max_coulomb_multipole_rank
         write (k) electron_count,csf_count,orbital_count
         write (k) diagonal_only,iccut(1),use_first_order_restriction
      end do

      do k = 30, 32+max_coulomb_multipole_rank
         close (k)
      end do

      return
      end subroutine setmcpmpi
