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
      program genmcpmpi
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 11 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 29 Jun 1998   *
!   Updated by Charlotte F. Fischer
!                                                                      *
!   Modified by Gediminas Gaigalas for new spin-angular integration.   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
      use parameter_def,   only:  max_orbitals
      use iso_fortran_env, only: int64
      use rangular_mpi_input, only: mpi_angular_input_options, parse_mpi_angular_options, &
           open_mpi_angular_profile, load_mpi_angular_header, load_mpi_angular_block, close_mpi_angular_input
      use rangular_profile, only: stage_clock, record_stage, close_stage_profile
      use memory_man
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use mpi_c
      use default_c, only: double, inpci, ndef, ndump
      use hblock_c, only: csf_counts_by_block, double, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use orb_c
      use stat_c, only: byte, csf_coupling_quantum_codes, csf_subshell_quantum_codes, double, jcupar, jqsar
      use mcp_c, only: diagonal_only, double, max_coulomb_multipole_rank, use_first_order_restriction
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i, only: getyn
      use setdbg_i, only: setdbg
      use setmc_i, only: setmc
      use setsum_i, only: setsum
      use cslh_i, only: cslh
      use strsum_i, only: strsum
      use factt_i, only: factt
      use settmp_i, only: settmp
      use mpiu_i, only: startmpi2
      use mpiu_i, only: stopmpi2
      use setdbgmpi_i, only: setdbgmpi
      use setmcpmpi_i, only: setmcpmpi
      use mcpmpi_i, only: mcpmpi
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nblk0 = 50   !-----------------------------------------------
      integer :: ncount1, ncore, nb
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical ::  debug, restrt, yes
      character(len=128) :: startdir, permdir, tmpdir
      character(len=8):: idblk(nblk0)
      character(len=3) :: idstring
      integer :: lenperm, lentmp
      type(mpi_angular_input_options) :: input_options
      integer(int64) :: program_start_count, stage_start_count
      character(len=16) :: input_format_name
!-----------------------------------------------
      call parse_mpi_angular_options(input_options)
      open(unit=739,file="rangular.log",status="UNKNOWN")
      startdir = " "
      permdir = " "
      tmpdir = " "
!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================
      call startmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length, ncount1, &
                           startdir, permdir, tmpdir, "RANGULAR_MPI")
      write (idstring, "(I3.3)") mpi_rank
      lenperm = len_trim (permdir)
      lentmp = len_trim (tmpdir)
      call open_mpi_angular_profile(input_options)
      input_format_name = "text"
      if (input_options%use_descriptors) input_format_name = "descriptors"
      program_start_count = stage_clock()
!=======================================================================
!  Get NDEF on node-0 and then send to all nodes
!=======================================================================
      if (mpi_rank == 0) then
         write (stderr_unit,*)
!         WRITE (istde,'(A)') ' Default settings?  (y/n) '
         write (stderr_unit,"(A)") " Full interaction?  (y/n) "
         yes = getyn ()
         if (yes) then
            ndef = 0
            write(739,"(A)") "y            ! Full interaction"
         else
            ndef = 1
            write(739,"(A)") "n            ! Full interaction"
         end if
      end if
      call mpi_barrier(mpi_comm_world,mpi_error_code)
      call mpi_bcast (ndef,1,mpi_integer,0,mpi_comm_world,mpi_error_code)

!=======================================================================
!
!  Checks and settings... Mostly done in backyard.
!
!    setdbg - debug output control parameters
!    setmc - machine- and installation- dependent constants
!    setsum - open the summary file
!    cslh - load header of the csl file
!    setmcp - open and check the  .mcp  files
!    strsum - append a summary of the inputs to the  .sum  file
!    factt - table of logarithms of factorials setup
!=======================================================================

      call setdbgmpi (debug, permdir(1:lenperm) // "/genmcp.dbg")
      call setmc
      stage_start_count = stage_clock()
      call load_mpi_angular_header(input_options,permdir,ncore,nblk0,idblk)
      call record_stage("input_header_and_broadcast",0,stage_start_count,storage=input_format_name)
      stage_start_count = stage_clock()
      restrt = .false.
!cjb  myid, nprocs = NOT args
!cjb  CALL setmcpmpi (myid, nprocs, ncore, idblk, 'mcp' // idstring)
      call setmcpmpi (ncore, idblk, "mcp" // idstring)
      if (ndef/=0 .and. mpi_rank==0) call strsum
      call factt
      call record_stage("mcp_setup_including_interaction_input",0,stage_start_count)

!=======================================================================
!     For each block, generate and sort the data
!=======================================================================

      do nb = 1, symmetry_block_count
         csf_count = csf_counts_by_block(nb)                        ! This ncf goes to common
         if (mpi_rank == 0) then
            write (6, *)
            write (6, *) "Block ", nb, ",  ncf = ", csf_count
         end if
            !*** Load current CSL block. Memories de-allocated in mcp ***
         stage_start_count = stage_clock()
         call alloc (csf_subshell_occupations, max_orbitals, csf_count, "IQA", "GENMCP")
         call alloc (csf_subshell_quantum_codes, max_orbitals,3,csf_count, "JQSA", "GENMCP")
         call alloc (csf_coupling_quantum_codes, max_orbitals, csf_count, "JCUPA", "GENMCP")
         call record_stage("csf_allocation",nb,stage_start_count,int(csf_count,int64))
!
         stage_start_count = stage_clock()
         call load_mpi_angular_block(input_options,ncore,nb)
         call record_stage("csf_load_and_broadcast",nb,stage_start_count,int(csf_count,int64),storage=input_format_name)
         !*** Open tmp.xx files for block nb ***
         call settmp (nb, max_coulomb_multipole_rank, "tmp" // idstring)
         !*** Generation of MCP coefficients ***
!cjb  myid, nprocs = NOT args
!cjb     CALL MCPmpi (NB, RESTRT, MYID, NPROCS, 'mcp' // idstring)
         call mcpmpi (nb, restrt, "mcp" // idstring)
      end do
      call close_mpi_angular_input(input_options)
      close(24)                                  ! Summary file
      close(739)                                 ! rangular.log
      if (debug) close(99)                       ! Debug file
      call record_stage("program_before_finalize_including_interaction_input",0,program_start_count)
      call close_stage_profile()
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================
      call stopmpi2 (mpi_rank,mpi_process_count,processor_name,processor_name_length,ncount1,"RANGULAR_MPI")
      stop
      end program genmcpmpi
