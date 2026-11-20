!***********************************************************************
      subroutine setmcp2(mpi_rank, mpi_process_count, core_subshell_count, block_labels, filehead)
      use mcp_generation_options, only: mcp_output_units
!
! A wrapper for setmcp/getinf. setmcp/getinf are then shared by serial
! and MPI programs.
!
!   Written by Xinghong He                Last revision: 30 Jun 1998
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use atomic_calculation_state, only: electron_count
      use csf_block_cutoffs, only: iccut, iccutblk
      use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
      use orbital_csf_basis, only: csf_count, orbital_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setmcp_i, only: setmcp
      use getinf_i, only: getinf
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: mpi_rank
      integer  :: mpi_process_count
      integer  :: core_subshell_count
      character  :: filehead*(*)
      character  :: block_labels(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k

      call setmcp (mpi_rank, mpi_process_count, core_subshell_count, block_labels, filehead)

! DIAG, ICCUTblk, LFORDR are set in GETINF

      if (mpi_rank == 0) call getinf

      do k = 30, 32 + max_coulomb_multipole_rank
         write (mcp_output_units(k)) electron_count, csf_count, orbital_count
         write (mcp_output_units(k)) diagonal_only, iccut(1), use_first_order_restriction
      end do

      do k = 30, 32 + max_coulomb_multipole_rank
         close(mcp_output_units(k))
      end do

      return
      end subroutine setmcp2
