
!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**        *****   ******  **   **  **    **   *****   ******        ***
!**       **   **  **      ***  **  ***  ***  **   **  **   **       ***
!**       **       **      ***  **  ** ** **  **       **   **       ***
!**       **  ***  ****    ** ****  ** ** **  **       ******        ***
!**       **   **  **      **  ***  **    **  **       **            ***
!**       **   **  **      **   **  **    **  **   **  **            ***
!**        *****   ******  **   **  **    **   *****   **            ***
!**                                                                  ***
!**   Program for generating the energy expression for H(DC). This   ***
!**   program is a derivative of GRASP2 (F. A. Parpia, I. P. Grant,  ***
!**   and C. F. Fischer, 1990).                                      ***
!**                                                                  ***
!**                         GRASP92 Version                          ***
!**          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
      program genmcp
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 11 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 29 Jun 1998   *
!   Updated by Charlotte F. Fischer                                    *
!                                                                      *
!   Modified by Gediminas Gaigalas:                         Feb 2017   *
!                               1) for new spin-angular integration,   *
!                               2) for sorting in the memory.          *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def,   only:  max_orbitals
      use iso_fortran_env, only: int64
      use rangular_input, only: angular_input_options, parse_angular_options, load_angular_header, &
                                load_angular_block, close_angular_input
      use rangular_profile, only: stage_clock, record_stage, close_stage_profile
      use memory_man
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use default_c
      use hblock_c
      use orb_c
      use stat_c
      use mcp_c
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use starttime_i
      use setdbg_i
      use setmc_i
      use setsum_i
      use cslh_i
      use setmcp2_i
      use strsum_i
      use factt_i
!      USE settmp_I
      use settmpgg_i
      use stoptime_i
      use lodcsh2_i, only: lodcsh2
      use mcp_i, only: mcp
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nblk0 = 50   !-----------------------------------------------
      integer :: ncount1, ncore, nb
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical ::  debug, restrt, yes
      character(len=128) :: startdir*128, permdir*128, tmpdir*128
      character(len=8):: idblk(nblk0)
      integer :: myid,nprocs
      type(angular_input_options) :: input_options
      integer(int64) :: program_start_count, stage_start_count
      character(len=16) :: input_format_name
!-----------------------------------------------------------------------
!
      call parse_angular_options(input_options)
      input_format_name = 'text'
      if (input_options%use_descriptors) input_format_name = 'descriptors'
      program_start_count = stage_clock()
      write(*,*)
      write(*,*) 'RANGULAR'
      write(*,*) 'This program performs angular integration '
      write(*,*) 'Input file: ', input_options%input_filename
      write(*,*) 'Outputfiles: mcp.30, mcp.31, ....'
      write(*,*) 'rangular.log                     '
      write(*,*)
      open(unit=739,file='rangular.log',status='UNKNOWN')
      myid = 0
      nprocs = 1
      call starttime (ncount1, 'RANGULAR')
!=======================================================================
!  Get NDEF
!=======================================================================
      if (myid == 0) then
!         WRITE (ISTDE, '(A)', ADVANCE='NO') 'Default settings?  (y/n) '
         write (stderr_unit,'(A)') ' Full interaction?  (y/n) '
         yes = getyn()
         if (yes) then
            ndef = 0
            write(739,'(A)') 'y            ! Full interaction'
         else
            ndef = 1
            write(739,'(A)') 'n            ! Full interaction'
         endif
      endif
!=======================================================================
!
!  Checks and settings... Mostly done in backyard.
!
!    chkplt - compatibility of plant substitutions
!    setdbg - debug output control parameters
!    setmc - machine- and installation- dependent constants
!    setsum - open the summary file
!    cslh - load header of the csl file
!    setmcp - open and check the  .mcp  files
!    strsum - append a summary of the inputs to the  .sum  file
!    factt - table of logarithms of factorials setup
!=======================================================================
!CFF  delete chkplt
!     CALL CHKPLT ('GENMCP')
      call setdbg (debug, 'genmcp.dbg')
      call setmc
!      IF (NDEF/=0 .AND. MYID==0) CALL SETSUM ('genmcp.sum')
      stage_start_count = stage_clock()
      call load_angular_header(input_options, ncore, nblk0, idblk)
      call record_stage('input_header', 0, stage_start_count, int(csf_count,int64), storage=input_format_name)
      restrt = .false.
      stage_start_count = stage_clock()
      call setmcp2 (myid, nprocs, ncore, idblk, 'mcp')
      if (ndef/=0 .and. myid==0) call strsum
      call factt
      call record_stage('setup_including_interaction_input', 0, stage_start_count)
!=======================================================================
!     For each block, generate and sort the data
!=======================================================================
      do nb = 1, symmetry_block_count
         csf_count = csf_counts_by_block(nb)                        ! This ncf goes to common
         if (myid == 0) then
            write (6, *)
            write (6, *) 'Block ', nb, ',  ncf = ', csf_count
         endif
            !*** Load current CSL block. Memories de-allocated in mcp ***
         stage_start_count = stage_clock()
         call alloc (csf_subshell_occupations, max_orbitals, csf_count, 'IQA', 'GENMCP')
         call alloc (csf_subshell_quantum_codes, max_orbitals,3,csf_count, 'JQSA', 'GENMCP')
         call alloc (csf_coupling_quantum_codes, max_orbitals, csf_count, 'JCUPA', 'GENMCP')
         call record_stage('csf_allocation', nb, stage_start_count)
!
         stage_start_count = stage_clock()
         call load_angular_block(input_options, ncore, nb)
         call record_stage('csf_load', nb, stage_start_count, int(csf_count,int64), storage=input_format_name)
         !*** Open tmp.xx files for block nb ***
!GG         CALL SETTMP (NB, KMAX, 'tmp')
         call settmpgg (nb, 30, 'tmp')
         !*** Generation of MCP coefficients ***
         call mcp (nb, restrt, myid, nprocs, 'mcp')
      end do
      call close_angular_input(input_options)
      call record_stage('program_total_including_interaction_input', 0, program_start_count)
      call close_stage_profile()
      close(24)                                  ! Summary file
      close(739)                                  ! rangular.log
      if (debug) close(99)                       ! Debug file
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================
      call stoptime (ncount1, 'RANGULAR')
      stop
      end program genmcp
