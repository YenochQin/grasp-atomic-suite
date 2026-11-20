!***********************************************************************
!                                                                      *
!cjb  myid, nprocs = NOT args
      subroutine setmcpmpi(ncore, idblk, filehead, descriptors)
!                                                                      *
! Set interaction options and write the MPI MCP file headers.          *
!                                                                      *
!   Written by Xinghong He                Last revision: 30 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!cjb mpi_C
!     USE mpi_state,    ONLY: MYID, NPROCS, ierr
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use csf_descriptor_reader, only: descriptor_source
      use mpi_csf_collectives, only: check_descriptor_status
      use csf_block_cutoffs, only: iccut
      use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
      use orbital_csf_basis, only: csf_count, orbital_count
      use atomic_calculation_state, only: electron_count
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
!cjb mpi_C
      use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_logical, mpi_rank
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
      type(descriptor_source), intent(in) :: descriptors
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!cjb mpi_C
!cjb  INTEGER :: K, ierr
      integer :: k, input_status
!-----------------------------------------------

      input_status = 0
      if (mpi_rank == 0) call getinf(descriptors, input_status)
      call check_descriptor_status(input_status)

      call setmcp (ncore, idblk, filehead)

      call mpi_bcast (iccut, size(iccut), mpi_integer, 0, mpi_comm_world, mpi_error_code)
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
