!***********************************************************************
!                                                                      *
      subroutine setcsl(name, ncore, idblk)
!                                                                      *
!  A container which calls lib92/lodcsh to get                         *
!     ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()                *
!  It then calls lodcsl to load the entire file (all blocks),          *
!  obtaining   IQA, JQSA, JCUPA                                        *
!  Note that unlike setcsl of mcp[mpi] and rci[mpi]vu, the following   *
!  quantities are assumed to be known:                                 *
!     nblock, ncfblk(), idblk(), ncftot.                               *
!  In this rscfmpi, they have been read into memory from mcp files via *
!  a call to setmcp.                                                   *
!                                                                      *
!  Xinghong He 98-08-06                                                *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def,    only: max_orbitals
      use memory_man, only: alloc, byte, double
      use hblock_c, only: csf_counts_by_block, double, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use def_c, only: double
      use orb_c, only: ncftot => csf_count, byte, csf_subshell_occupations, double, max_orbitals
      use stat_c, only: byte, csf_coupling_quantum_codes, csf_subshell_quantum_codes, double, jcupar, jqsar
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i, only: openfl
      use lodcsh_i, only: lodcsh
      use lodcsh2gg_i, only: lodcsh2gg
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ncore
      character  :: name*(*)
      character  :: idblk(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer ::  ierr, ios
      character :: str*15
!-----------------------------------------------
!
! opens, reads the header of the file

      call openfl (21, name, "FORMATTED", "OLD", ierr)
      if (ierr == 1) then
         write (6, *) "Error when opening ", name(1:len_trim(name))
         error stop
      end if

      read (21, "(1A15)", iostat=ios) str
      if (ios/=0 .or. str/="Core subshells:") then
         write (6, *) "Not a Configuration Symmetry List File;"
         close(21)
         error stop
      end if

      !..Load header of <name> file
      call lodcsh (21, ncore)

! Allocate memories for all blocks and then load the entire file

      call alloc (csf_subshell_occupations, max_orbitals, ncftot, "IQA", "SETCSL")
!GG      CALL ALLOC (jqsa, NNNW,3,NCFTOT, 'JQSA', 'SETCSL')
!GG      CALL ALLOC (jcupa, NNNW, NCFTOT, 'JCUPA', 'SETCSL')

      call lodcsh2gg (21, ncore, -119)
      ! -119 means "load all blocks"

      return
      end subroutine setcsl
