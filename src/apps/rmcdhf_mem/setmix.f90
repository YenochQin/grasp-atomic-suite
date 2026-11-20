!***********************************************************************
      subroutine setmix(name)
!
!   Opens the  .mix  file on stream 25; writes a header to this file.  *
!                                                                      *
!   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
!   Modified by Xinghong He               Last revision: 13 Jul 1998   *
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use atomic_calculation_state, only: electron_count
      use first_order_csf_cutoffs, only: iccut
      use mcp_file_metadata, only: max_mcp_file_multipole_rank
      use rmcdhf_mcp_options, only: diag, lfordr
      use orbital_csf_basis, only: csf_count, orbital_count
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use file_open, only: openfl
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character  :: name*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: open_status
!-----------------------------------------------
!     POINTER (PCCMIN,ICCMIN(1))
!     POINTER (PNTRIQ,RIQDUM)
!
!
!     !...Wants nblock only
!     POINTER (pncfblk,ncfblk(0:*))
!
!-----------------------------------------------------------------------
      call openfl (25, name, "UNFORMATTED", "NEW", open_status)
      if (open_status /= 0) then
         write (stderr_unit, *) "Error when opening ", name(1:len_trim(name))
         stop
      end if
!
!   Write the file header
!
      write (25) "G92MIX"
      write (25) electron_count, csf_count, orbital_count, 0, 0, symmetry_block_count
!     ...The above record will be overidden in matrix.f
!        with the final form of
!     WRITE (25)  NELEC, NCF, NW, nvectot, nvecsiz, nblock

      return
      end subroutine setmix
