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
      program rscfvu
!                                                                      *
!   Entry routine for RSCFVU. Controls the entire computation.         *
!                                                                      *
!   Call(s) to: [LIB92]: SETMC, SETCON.                                *
!               [RSCF92]: CHKPLT, setcsl, setdbg                       *
!                        SETMIX, SETRES, SETSUM, STRSUM.               *
!               [NJGRAF]: FACTT.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 31 Dec 1992   *
!   Block version by Xinghong He          Last revision: 17 Aug 1998   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
       use calculation_input_defaults, only: inpci, ndef, ndump
       use closed_core_subshells, only: closed_core_subshell_count
       use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
       use mpi_s, only: mpi_process_count, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yes_no_input, only: getyn
      use run_timer_start, only: starttime
      use setdbg_i, only: setdbg
      use machine_constants, only: setmc
      use physical_constants, only: setcon
      use setsum_i, only: setsum
      use setmcp_i, only: setmcp
      use setcsl_i, only: setcsl
      use getscd_i, only: getscd
      use strsum_i, only: strsum
      use setmix_i, only: setmix
      use log_factorial_table, only: factt
      use scf_i, only: scf
      use run_timer_finish, only: stoptime
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nblk0 = 50  ! Maximum number of blocks
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncore1, ncount1
      logical :: eol, yes
      character :: idblk(nblk0)*8
!-----------------------------------------------
!
!
! Things for timing
!-----------------------------------------------------------------------
      mpi_rank = 0
      mpi_process_count = 1

      write(*,*)
      write(*,*) "RMCDHF"
      write(*,*) "This program determines the radial orbitals   "
      write(*,*) "and the expansion coefficients of the CSFs         "
      write(*,*) "in a self-onsistent field proceedure               "
      write(*,*) "Input file:  isodata, rcsf_descriptors.parquet + rcsf_header.toml, rwfn.inp, mcp.30, ..."
      write(*,*)                                                       &
               "Outputfiles: rwfn.out, rmix.out, rmcdhf.sum, rmcdhf.log"
      write(*,*)

      call starttime (ncount1, "RMCDHF")

     open(unit=734,file="rmcdhf.log",status="UNKNOWN")

!=======================================================================
!  Get NDEF
!=======================================================================

      if (mpi_rank == 0) then
         write (stderr_unit, "(A)", advance="NO") "Default settings?  (y/n) "
         yes = getyn()
         if (yes) then
            ndef = 0
            write(734,"(A)") "y            ! Default settings"
         else
            ndef = 1
         end if
      end if

!=======================================================================
!
!  Checks and settings... Mostly done in backyard.
!
!    SETDBG - Debug output control parameters
!    SETMC - machine- and installation- dependent constants
!    SETCON - physical constants
!    SETSUM - open the summary file
!    SETMCP - open and check the  .mcp  files
!    SETCSL - open, check, load data from, and close the  .csl  file
!    STRSUM - append a summary of the inputs to the  .sum  file
!    SETMIX - mixing coefficients file setup
!    FACTT - table of logarithms of factorials setup
!=======================================================================

      call setdbg ("rscf92.dbg")
      call setmc
      call setcon

      call setsum ("rmcdhf.sum")

      call setmcp (closed_core_subshell_count, nblk0, idblk, "mcp")
      call setcsl ("rcsf_descriptors.parquet", ncore1, idblk)
      if (closed_core_subshell_count /= ncore1) stop "rscfvu: ncore"

!=======================================================================
!  Gather all remaining information and perform some setup. This
!  part (routine) asks for user-inputs.
!=======================================================================

      call getscd (eol, idblk, "isodata", "rwfn.inp", "npol.inp")

      if (mpi_rank == 0) then
         call strsum
         if (eol) call setmix ("rmix.out")
      end if

      call factt

!=======================================================================
!  Proceed with the SCF calculation close all files except
!  the  .sum  file
!=======================================================================

      call scf (eol, "rwfn.out")
      close (734)

!=======================================================================
!  Execution finished; Statistics output
!=======================================================================

      call stoptime (ncount1, "RMCDHF")

      stop
      end program rscfvu
