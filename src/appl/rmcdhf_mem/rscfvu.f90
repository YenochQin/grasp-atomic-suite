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
       use default_c
       use core_c
       use iounit_c
       use mpi_s
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use starttime_i
      use setdbg_i
      use setmc_i
      use setcon_i
      use setsum_i
      use setmcp_i
      use setcsl_i
      use getscd_i
      use strsum_i
      use setmix_i
      use factt_i
      use scf_i
      use stoptime_i
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
      write(*,*) 'RMCDHF'
      write(*,*) 'This program determines the radial orbitals   '
      write(*,*) 'and the expansion coefficients of the CSFs         '
      write(*,*) 'in a self-onsistent field proceedure               '
      write(*,*) 'Input file:  isodata, rcsf.inp, rwfn.inp, mcp.30, ...'
      write(*,*)                                                       &
               'Outputfiles: rwfn.out, rmix.out, rmcdhf.sum, rmcdhf.log'
      write(*,*)

      call starttime (ncount1, 'RMCDHF')

     open(unit=734,file='rmcdhf.log',status='UNKNOWN')

!=======================================================================
!  Get NDEF
!=======================================================================

      if (mpi_rank == 0) then
         write (stderr_unit, '(A)', advance='NO') 'Default settings?  (y/n) '
         yes = getyn()
         if (yes) then
            ndef = 0
            write(734,'(A)') 'y            ! Default settings'
         else
            ndef = 1
         endif
      endif

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

      call setdbg ('rscf92.dbg')
      call setmc
      call setcon

      call setsum ('rmcdhf.sum')

      call setmcp (closed_core_subshell_count, nblk0, idblk, 'mcp')
      call setcsl ('rcsf.inp', ncore1, idblk)
      if (closed_core_subshell_count /= ncore1) stop 'rscfvu: ncore'

!=======================================================================
!  Gather all remaining information and perform some setup. This
!  part (routine) asks for user-inputs.
!=======================================================================

      call getscd (eol, idblk, 'isodata', 'rwfn.inp', 'npol.inp')

      if (mpi_rank == 0) then
         call strsum
         if (eol) call setmix ('rmix.out')
      endif

      call factt

!=======================================================================
!  Proceed with the SCF calculation close all files except
!  the  .sum  file
!=======================================================================

      call scf (eol, 'rwfn.out')
      close (734)

!=======================================================================
!  Execution finished; Statistics output
!=======================================================================

      call stoptime (ncount1, 'RMCDHF')

      stop
      end program rscfvu
