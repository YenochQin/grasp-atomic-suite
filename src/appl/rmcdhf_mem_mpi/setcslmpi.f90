!***********************************************************************
!
      subroutine setcslmpi(name, ncore, idblk)
!
!  A container which calls lib92/lodcsh to get
!     ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
!  It then calls lodcsl to load the entire file (all blocks), obtaining
!     IQA, JQSA, JCUPA
!  Note that unlike setcsl of mcp[mpi] and rci[mpi]vu, the following
!  quantities are assumed to be known:
!     nblock, ncfblk(), idblk(), ncftot.
!  In this rscfmpi, they have been read into memory from mcp files via
!  a call to setmcp.
!
!  Xinghong He 98-08-06
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNWP*NCF),                  *
!                                  JCUPA(NNNWP*NCF)                    *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def,    only: max_orbitals
      use memory_man
      use hblock_c
      use def_c
      use orb_c, ncftot=>csf_count
      use stat_c
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      use lodcslmpigg_i
      use lodcsh_i, only: lodcsh
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
      integer ::  ios
      character :: str*15
!-----------------------------------------------
!
! node-0 opens, reads the header of the file

      if (mpi_rank .eq. 0) then
         call openfl (21, name, 'FORMATTED', 'OLD', mpi_error_code)
         if (mpi_error_code .eq. 1) then
            print *, 'Error when opening ',name(1:len_trim (name))
            stop
         endif

         read (21,'(1A15)',iostat = ios) str
         if ((ios .ne. 0) .or.                   &
               (str .ne. 'Core subshells:')) then
            print *, 'Not a Configuration Symmetry List File;'
            close (21)
            stop
         endif

         !..Load header of <name> file
         call lodcsh (21, ncore)
      end if

! Broadcast results to other nodes. ncfblk should be allocated
! on these nodes (with myid .ne. 0)

      call mpi_bcast (orbital_count,    1, mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (ncore, 1, mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (electron_count, 1, mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (orbital_principal_n(1),orbital_count, mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (orbital_kappa(1),orbital_count, mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (orbital_l(1),orbital_count, mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (orbital_two_j(1),orbital_count, mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (orbital_symmetry_labels(1),2*orbital_count, mpi_character,0,mpi_comm_world,mpi_error_code)

! Allocate memories for all blocks and then load the entire file

      call alloc (csf_subshell_occupations, max_orbitals, ncftot, 'IQA', 'SETCSLmpi')
!GG      CALL ALLOC (jqsa, NNNW,3,NCFTOT, 'JQSA', 'SETCSmpiL')
!GG      CALL ALLOC (jcupa, NNNW, NCFTOT, 'JCUPA', 'SETCSLmpi')

!GG      CALL LODCSLmpi (21, NCORE, -119)
      call lodcslmpigg (21, ncore, -119)
      ! -119 means "load all blocks"

      return
      end subroutine setcslmpi
