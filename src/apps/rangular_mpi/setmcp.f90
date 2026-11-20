!***********************************************************************
!                                                                      *
!cjb  myid, nprocs = NOT args
      subroutine setmcp(core_subshell_count, block_labels, filehead)
!                                                                      *
!   Open and check the  .mcp  files. File 30 stores the structure of   *
!   H(DC) ; file 31 stores the  T  coefficients;  files 32, 33, ...,   *
!   store V(0), V(1), ... .                                            *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT, GETYN, OPENFL.                        *
!               [GENMCP]: GETINF.                                      *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 08 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 30 Jun 1998   *
!
!   Used by mcpvu, mcpmpivu
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use mcp_generation_options,     only: max_coulomb_multipole_rank
      use orbital_csf_basis,     only: orbital_count, orbital_two_j
      use hamiltonian_symmetry_blocks,  only: symmetry_block_count, csf_counts_by_block
      use standard_io_units,  only: stderr_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use integer_format, only: convrt
      use file_open, only: openfl
      use mpi_state,     only: mpi_rank, mpi_process_count, mpi_error_code
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!     INTEGER, INTENT(IN) :: MYID
!     INTEGER, INTENT(IN) :: NPROCS
      integer, intent(in) :: core_subshell_count
      character, intent(in) :: filehead*(*)
      character, intent(in) :: block_labels(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, lng, text_length, i
      logical :: found, found1, getyn, yes
      character :: number_text*2
!-----------------------------------------------
!   Determine KMAX; this is the number of  .mcp  files for the
!   two-electron integrals

      max_coulomb_multipole_rank = 0
      do k = 1, orbital_count
         max_coulomb_multipole_rank = max(max_coulomb_multipole_rank,orbital_two_j(k))
      end do

!   All files  mcp.xx  are UNFORMATTED;

      lng = len_trim(filehead)
      do k = 30, 32 + max_coulomb_multipole_rank
         call convrt (k, number_text, text_length)
         call openfl (k, filehead(1:lng)//"."//number_text(1:2), "UNFORMATTED", &
            "UNKNOWN", mpi_error_code)
         if (mpi_error_code == 0) cycle
         do i = 30, k
            close(i)
         end do
         write (stderr_unit, *) "Error when opening the mcp files"
         stop
      end do
!
!   We want to know kmax before openning other mcp files (not mcp.30)
!   in rscf
!
      write (30) core_subshell_count, symmetry_block_count, max_coulomb_multipole_rank
      write (30) (csf_counts_by_block(i),i=1,symmetry_block_count)
      write (30) (block_labels(i),i=1,symmetry_block_count)

      do k = 30, 32 + max_coulomb_multipole_rank
         write (k) "MCP", symmetry_block_count, mpi_rank, mpi_process_count
      end do

      return
      end subroutine setmcp
