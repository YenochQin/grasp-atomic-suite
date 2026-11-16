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
      use def_c
      use foparm_c, only: iccut
      use mcpa_c, only: double, max_mcp_file_multipole_rank
      use mcpb_c, only: diag, double, lfordr
      use orb_c
      use hblock_c, only: csf_counts_by_block, double, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i, only: openfl
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character  :: name*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ierr
!-----------------------------------------------
!     POINTER (PCCMIN,ICCMIN(1))
!     POINTER (PNTRIQ,RIQDUM)
!
!
!     !...Wants nblock only
!     POINTER (pncfblk,ncfblk(0:*))
!
!-----------------------------------------------------------------------
      call openfl (25, name, "UNFORMATTED", "NEW", ierr)
      if (ierr /= 0) then
         write (stderr_unit, *) "Error when opening ", name(1:len_trim(name))
         error stop
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
