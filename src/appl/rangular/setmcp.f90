!***********************************************************************
!                                                                      *
      subroutine setmcp(myid, nprocs, ncore, idblk, filehead)
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
!   M o d u l e s
!-----------------------------------------------
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
      use orbital_csf_basis, only: orbital_count, orbital_two_j
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use integer_format, only: convrt
      use file_open, only: openfl
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: myid
      integer, intent(in) :: nprocs
      integer, intent(in) :: ncore
      character, intent(in) :: filehead*(*)
      character, intent(in) :: idblk(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, lng, lck, ierr, i
      logical :: found, found1, getyn, yes
      character :: ck*2
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
         call convrt (k, ck, lck)
         call openfl (k, filehead(1:lng)//"."//ck(1:2), "UNFORMATTED", &
            "UNKNOWN", ierr)
         if (ierr == 0) cycle
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
      write (30) ncore, symmetry_block_count, max_coulomb_multipole_rank
      write (30) (csf_counts_by_block(i),i=1,symmetry_block_count)
      write (30) (idblk(i),i=1,symmetry_block_count)

      do k = 30, 32 + max_coulomb_multipole_rank
         write (k) "MCP", symmetry_block_count, myid, nprocs
      end do

      return
      end subroutine setmcp
