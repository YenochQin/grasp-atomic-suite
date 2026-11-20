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
      use, intrinsic :: iso_fortran_env, only: output_unit
      use calculation_file_units, only: angular_log_unit, summary_output_unit
      use debug_output_settings, only: angular_debug_unit
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 11 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 29 Jun 1998   *
!   Updated by Charlotte F. Fischer
!                                                                      *
!   Modified by Gediminas Gaigalas for new spin-angular integration.   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
      use mpi_csf_collectives, only: allocate_descriptor_arrays
      use grasp_parameters,   only:  max_orbitals
      use iso_fortran_env, only: int64
      use rangular_mpi_input, only: mpi_angular_input_options, parse_mpi_angular_options, &
           open_mpi_angular_profile, load_mpi_angular_header, load_mpi_angular_block, close_mpi_angular_input
      use rangular_profile, only: stage_clock, record_stage, close_stage_profile
      use array_memory_management, only: alloc
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use mpi_state, only: mpi_process_count, mpi_rank, processor_name, processor_name_length
      use calculation_input_defaults, only: inpci, ndef, ndump
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use orbital_csf_basis, only: csf_count, csf_subshell_occupations, max_orbitals
      use csf_coupling_data, only: csf_coupling_quantum_codes, csf_subshell_quantum_codes, jcupar, jqsar
      use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use calculation_state_files, only: state_files, state_name_length, files_for_state
      use mpi_calculation_state, only: select_mpi_state
      use setdbg_i, only: setdbg
      use machine_constants, only: setmc
      use setsum_i, only: setsum

      use strsum_i, only: strsum
      use log_factorial_table, only: factt
      use settmp_i, only: settmp
      use mpi_utilities, only: startmpi2
      use mpi_utilities, only: stopmpi2
      use setdbgmpi_i, only: setdbgmpi
      use setmcpmpi_i, only: setmcpmpi
      use mcpmpi_i, only: mcpmpi
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nblk0 = 50   !-----------------------------------------------
      integer :: start_clock_count, core_subshell_count, symmetry_block_index
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical ::  debug, is_restart, yes
      character(len=128) :: startup_directory, permdir, tmpdir
      character(len=8):: idblk(nblk0)
      character(len=3) :: block_label_text
      type(mpi_angular_input_options) :: input_options
      type(state_files) :: files
      character(state_name_length) :: state
      integer(int64) :: program_start_count, stage_start_count
      character(len=16) :: input_format_name
!-----------------------------------------------
      call parse_mpi_angular_options(input_options)
      startup_directory = " "
      permdir = " "
      tmpdir = " "
!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================
      call startmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length, start_clock_count, &
                           startup_directory, permdir, tmpdir, "RANGULAR_MPI")
      write (block_label_text, "(I3.3)") mpi_rank
      call select_mpi_state(state, "Use zero-first partition from descriptor header? (y/n)", yes, permdir)
      files = files_for_state(state)
      if (input_options%input_filename == "") input_options%input_filename = files%descriptors
      if (mpi_rank == 0) then
         open(newunit=angular_log_unit,file=trim(startup_directory)//"/"//files%angular_log,status="UNKNOWN")
         if (state /= "") write(angular_log_unit,"(A)") trim(state)
      end if
      call open_mpi_angular_profile(input_options)
      input_format_name = "descriptors"
      program_start_count = stage_clock()
!=======================================================================
!  Get NDEF on node-0 and then send to all nodes
!=======================================================================
      ndef = merge(1, 0, yes)
      if (mpi_rank == 0) write(angular_log_unit,"(A)") &
         merge("y", "n", yes)//"            ! Use zero-first partition from descriptor header"

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

      call setdbgmpi (debug, trim(permdir)//"/"//files%angular_debug)
      call setmc
      stage_start_count = stage_clock()
      call load_mpi_angular_header(input_options,permdir,core_subshell_count,nblk0,idblk)
      call record_stage("input_header_and_broadcast",0,stage_start_count,storage=input_format_name)
      stage_start_count = stage_clock()
      is_restart = .false.
!cjb  myid, nprocs = NOT args
!cjb  CALL setmcpmpi (myid, nprocs, ncore, idblk, 'mcp' // idstring)
      call setmcpmpi (core_subshell_count, idblk, files%mcp, input_options%descriptors)
      if (ndef/=0 .and. mpi_rank==0) call strsum
      call factt
      call record_stage("mcp_setup_including_interaction_input",0,stage_start_count)

!=======================================================================
!     For each block, generate and sort the data
!=======================================================================

      do symmetry_block_index = 1, symmetry_block_count
         csf_count = csf_counts_by_block(symmetry_block_index)                        ! This ncf goes to common
         if (mpi_rank == 0) then
            write (output_unit, *)
            write (output_unit, *) "Block ", symmetry_block_index, ",  ncf = ", csf_count
         end if
            !*** Load current CSL block. Memories de-allocated in mcp ***
         stage_start_count = stage_clock()
         call allocate_descriptor_arrays(csf_count,csf_subshell_occupations, &
                                         csf_subshell_quantum_codes,csf_coupling_quantum_codes)
         call record_stage("csf_allocation",symmetry_block_index,stage_start_count,int(csf_count,int64))
!
         stage_start_count = stage_clock()
         call load_mpi_angular_block(input_options,core_subshell_count,symmetry_block_index)
         call record_stage("csf_load_and_broadcast",symmetry_block_index,stage_start_count,int(csf_count,int64), &
            & storage=input_format_name)
         !*** Open tmp.xx files for block nb ***
         call settmp (symmetry_block_index, max_coulomb_multipole_rank, files%temporary//block_label_text)
         !*** Generation of MCP coefficients ***
!cjb  myid, nprocs = NOT args
!cjb     CALL MCPmpi (NB, RESTRT, MYID, NPROCS, 'mcp' // idstring)
         call mcpmpi (symmetry_block_index, is_restart, files%mcp)
      end do
      call close_mpi_angular_input(input_options)
      if (summary_output_unit /= -1) close(summary_output_unit)                                  ! Summary file
      if (mpi_rank == 0) close(angular_log_unit)                                 ! rangular.log
      if (debug .and. angular_debug_unit /= -1) close(angular_debug_unit)                       ! Debug file
      call record_stage("program_before_finalize_including_interaction_input",0,program_start_count)
      call close_stage_profile()
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================
      call stopmpi2 (mpi_rank,mpi_process_count,processor_name,processor_name_length,start_clock_count,"RANGULAR_MPI")
      stop
      end program genmcpmpi
