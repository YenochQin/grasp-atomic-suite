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
       use default_c
       use core_c
       use iounit_c
       use mpi_c
!GG po to isimti
       use tatb_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use setdbgmpi_i
      use setmc_i
      use setcon_i
      use setsum_i
      use setmcp_i
      use setcslmpi_i
      use getscdmpi_i
      use strsum_i
      use setmix_i
      use factt_i
      use scfmpi_i
      use npol_mpi, only: npol_abort_mpi
      use mpiu_i, only: startmpi2
      use mpiu_i, only: stopmpi2
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nblk0 = 50  ! Maximum number of blocks
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncore1, ncount1, lenperm, lentmp
      logical :: eol, yes
      character :: idblk(nblk0)*8
!      CHARACTER (LEN = *) :: host
      character(:), allocatable :: npfile
      character (len = 3) :: idstring
      character (len = 128) :: startdir,permdir,tmpdir,file_rcsl,file1,file2
!-----------------------------------------------
!
      open(unit=734,file='rmcdhf.log',status='UNKNOWN')

!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================
      startdir = '  '  ;    file_rcsl = '  '
      permdir = '  '   ;    file1     = '  '
      tmpdir = '  '    ;    file2     = '  '
      call startmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length, ncount1,       &
                      startdir, permdir, tmpdir, 'RMCDHF_MPI')
      write (idstring, '(I3.3)') mpi_rank
      lenperm = len_trim (permdir)
      lentmp = len_trim (tmpdir)

!=======================================================================
!  Get NDEF
!=======================================================================

      if (mpi_rank == 0) then
         write (stderr_unit,*)
         write (stderr_unit, '(A)') 'Default settings?  (y/n) '
         yes = getyn()
         if (yes) then
            ndef = 0
!cjb fort.734
!cjb        WRITE(734,'(A)') 'y            ! Default settings'
         else
            ndef = 1
         endif
      endif
     call mpi_bcast (ndef,1,mpi_integer,0,mpi_comm_world,mpi_error_code)

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
      call setdbgmpi (permdir(1:lenperm) // '/rscf92.dbg')
      call setmc
      call setcon

      if (mpi_rank .eq. 0) call setsum (permdir(1:lenperm) // '/rmcdhf.sum')

      call setmcp (closed_core_subshell_count, nblk0, idblk, 'mcp' // idstring)
      if(mpi_rank == 0) file_rcsl = permdir(1:lenperm)//'/rcsf.inp'
      call setcslmpi (file_rcsl, ncore1, idblk)
      if (closed_core_subshell_count /= ncore1) stop 'rscfmpivu: ncore'

!=======================================================================
!  Gather all remaining information and perform some setup. This
!  part (routine) asks for user-inputs.
!=======================================================================

      if(mpi_rank == 0) file1 = permdir(1:lenperm) // '/isodata'
      if(mpi_rank == 0) file2 = permdir(1:lenperm) // '/rwfn.inp'
      npfile = ''
      if (mpi_rank == 0) then
         if (lenperm == len(permdir)) call npol_abort_mpi(mpi_comm_world, &
              'npol: permanent input directory reaches the legacy path limit')
         npfile = trim(permdir) // '/npol.inp'
      end if
      call getscdmpi (eol, idblk, file1, file2, npfile)
!     &            permdir(1:lenperm) // '/isodata',
!     &            permdir(1:lenperm) // '/rwfn.inp')
      file1 = '  '
      if(mpi_rank == 0) file1 = permdir(1:lenperm) // '/rmix.out'
      if (mpi_rank .eq. 0) then
         call strsum
         if (eol) call setmix (file1)
      endif

      call factt

!=======================================================================
!  Proceed with the SCF calculation close all files except
!  the  .sum  file
!=======================================================================

      file1 = '  '
      if(mpi_rank == 0) file1 = permdir(1:lenperm) // '/rwfn.out'

      call scfmpi (eol, file1)

!=======================================================================
!  Execution finished; Statistics output
!=======================================================================

      call stopmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length, ncount1, 'RMCDHF_MPI')

      end program rscfmpivu
