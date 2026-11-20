!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**       ******    *****    *****   *******   *****    *****        ***
!**       **   **  **   **  **   **  **       **   **  **   **       ***
!**       **   **  **       **       **       **   **       **       ***
!**       ******    *****   **       ****      *****       **        ***
!**       **  **        **  **       **          **       **         ***
!**       **   **  **   **  **   **  **         **      **           ***
!**       **   **   *****    *****   **        **      *******       ***
!**                                                                  ***
!**            Relativistic Self-Consistent-Field Program            ***
!**                                                                  ***
!**   This program is a derivative of GRASP2 (F. A. Parpia, I. P.    ***
!**   Grant, and C. F. Fischer, 1990)                                ***
!**                                                                  ***
!**                            GRASP92                               ***
!**          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
!**   Modified by C.F. FIscher for block input                       ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
      program rscfmpivu
!                                                                      *
!   Entry routine for RSCFVU. Controls the entire computation.      *
!                                                                      *
!   Call(s) to: [LIB92]: SETMC, SETCON.                                *
!               [RSCF92]: CHKPLT, setcsl, setdbg
!                        SETMIX, SETRES, SETSUM, STRSUM.
!               [NJGRAF]: FACTT.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 31 Dec 1992   *
!   Block version by Xinghong He          Last revision: 17 Aug 1998   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNWP*NCF),                  *
!                                  JCUPA(NNNWP*NCF)                    *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
       use calculation_input_defaults, only: inpci, ndef, ndump
       use closed_core_subshells, only: closed_core_subshell_count
       use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
       use mpi_state, only: mpi_comm_world, mpi_process_count, mpi_rank, &
          & processor_name, processor_name_length
!GG po to isimti
       use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use calculation_state_files, only: check_orbital_output
      use calculation_state_files, only: state_files, state_name_length, files_for_state
      use mpi_calculation_state, only: select_mpi_state
      use setdbgmpi_i, only: setdbgmpi
      use machine_constants, only: setmc
      use physical_constants, only: setcon
      use setsum_i, only: setsum
      use setmcp_i, only: setmcp
      use setcslmpi_i, only: setcslmpi
      use getscdmpi_i, only: getscdmpi
      use strsum_i, only: strsum
      use setmix_i, only: setmix
      use log_factorial_table, only: factt
      use scfmpi_i, only: scfmpi
      use mpi_nuclear_polarization, only: npol_abort_mpi
      use mpi_utilities, only: startmpi2
      use mpi_utilities, only: stopmpi2
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nblk0 = 50  ! Maximum number of blocks
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncore1, ncount1, lenperm
      type(state_files) :: files
      character(state_name_length) :: state
      logical :: eol, yes
      character :: idblk(nblk0)*8
!      CHARACTER (LEN = *) :: host
      character(:), allocatable :: npfile
      character (len = 128) :: startdir,permdir,tmpdir
      character(:), allocatable :: file_rcsl, file1, file2
!-----------------------------------------------
!

!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================
      startdir = "  "  ;    file_rcsl = "  "
      permdir = "  "   ;    file1     = "  "
      tmpdir = "  "    ;    file2     = "  "
      call startmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length, ncount1,       &
                      startdir, permdir, tmpdir, "RMCDHF_MPI")
      lenperm = len_trim (permdir)

!=======================================================================
!  Get NDEF
!=======================================================================

      call select_mpi_state(state, "Default settings? (y/n)", yes, permdir)
      files = files_for_state(state)
      if (mpi_rank == 0) call check_orbital_output(trim(permdir)//"/"//files%orbitals, writing_initial=.false.)
      ndef = merge(0, 1, yes)
      if (mpi_rank == 0) then
         open(unit=734,file=trim(startdir)//"/"//files%scf_log,status="UNKNOWN")
         if (state /= "") then
            write(734,"(A)") trim(state)
            write(734,"(A)") merge("y","n",yes)//"            ! Default settings"
         end if
      end if

!=======================================================================
!
!  Checks and settings... Mostly done in backyard.
!
!    CHKPLT - compatibility of plant substitutions
!    SETDBGmpi - Debug output control parameters
!    SETMC - machine- and installation- dependent constants
!    SETCON - physical constants
!    SETSUM - open the summary file
!    SETMCP - open and check the  .mcp  files
!    SETCSLmpi - open, check, load data from, and close the  .csl  file
!    STRSUM - append a summary of the inputs to the  .sum  file
!    SETMIX - mixing coefficients file setup
!    FACTT - table of logarithms of factorials setup
!=======================================================================
      call setdbgmpi (trim(permdir)//"/"//files%scf_debug)
      call setmc
      call setcon

      if (mpi_rank == 0) call setsum (trim(permdir)//"/"//files%scf_summary)

      call setmcp (closed_core_subshell_count, nblk0, idblk, files%mcp)
      if(mpi_rank == 0) file_rcsl = trim(permdir)//"/"//files%descriptors
      call setcslmpi (file_rcsl, ncore1, idblk)
      if (closed_core_subshell_count /= ncore1) stop "rscfmpivu: ncore"

!=======================================================================
!  Gather all remaining information and perform some setup. This
!  part (routine) asks for user-inputs.
!=======================================================================

      if(mpi_rank == 0) file1 = permdir(1:lenperm) // "/isodata"
      if(mpi_rank == 0) file2 = trim(permdir)//"/"//files%initial_orbitals
      npfile = ""
      if (mpi_rank == 0) then
         if (lenperm == len(permdir)) call npol_abort_mpi(mpi_comm_world, &
              "npol: permanent input directory reaches the legacy path limit")
         npfile = trim(permdir) // "/npol.inp"
      end if
      call getscdmpi (eol, idblk, file1, file2, npfile)
      file1 = "  "
      if(mpi_rank == 0) file1 = trim(permdir)//"/"//files%mixing
      if (mpi_rank == 0) then
         call strsum
         if (eol) call setmix (file1)
      end if

      call factt

!=======================================================================
!  Proceed with the SCF calculation close all files except
!  the  .sum  file
!=======================================================================

      file1 = "  "
      if(mpi_rank == 0) file1 = trim(permdir)//"/"//files%orbitals

      call scfmpi (eol, file1)

!=======================================================================
!  Execution finished; Statistics output
!=======================================================================

      call stopmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length, ncount1, "RMCDHF_MPI")

      end program rscfmpivu
