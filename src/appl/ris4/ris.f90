!***********************************************************************
!                                                                      *
      program ris
!                                                                      *
!   Entry routine for RIS. Controls the entire computation.       *
!                                                                      *
!   Call(s) to: [LIB92]: GETMIX, SETCSL, SETMC, SETCON.                *
!               [SMS92]: CHKPLT, GETSMD, SETDBG, SETSUM                *
!                        STRSUM                                        *
!               [NJGRAF]: FACTT.                                       *
!                                                                      *
!   Written by Per Jonsson                Last revision: 17 Jan 1996   *
!   Modify  by Gediminas Gaigalas                        26 Oct 2009   *
!   Modified by J. Ekman                                 25 Nov 2013   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use calculation_input_defaults, only: inpci, ndef, ndump
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use debug_output_settings,         only: cutoff
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use calculation_state_files, only: read_state_name
      use yes_no_input, only: getyn
      use setdbg_i, only: setdbg
      use machine_constants, only: setmc
      use physical_constants, only: setcon
      use setsum_i, only: setsum
      use csf_basis_loader, only: setcsla
      use getmixblock_i, only: getmixblock
      use strsum_i, only: strsum
      use log_factorial_table, only: factt
      use ris_cal_i, only: ris_cal
      use getsmd_i, only: getsmd
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: dr2
      logical   :: yes
      character :: name*24
      integer :: nci, ncore_not_used, nopar
!-----------------------------------------------
!
!   Matrix elements smaller than CUTOFF are not accumulated
!
      cutoff = 1.0e-10_real64
!
      write (stderr_unit, *)
      write (stderr_unit, *) "RIS: Execution begins ..."
!      CALL STARTTIME (ncount1, 'RIS')

      write (stderr_unit, *)
      write (stderr_unit, *) "Default settings?"
      yes = getyn ()
      write (stderr_unit, *)
      if (yes) then
         ndef = 0
      else
         ndef = 1
      end if

      call read_state_name(name)
      print *

      write (stderr_unit, *) "Mixing coefficients from a CI calc.?"
      yes = getyn ()
      if (yes) then
         nci = 0
      else
         nci = 1
      end if
      write (stderr_unit, *)

!
!   Check compatibility of plant substitutions
!
!GG      CALL CHKPLT
!
!   Determine if there is to be any debug printout; this will be
!   made on the  .dbg  file
!
      call setdbg
!
!   Perform machine- and installation-dependent setup
!
      call setmc
!
!   Set up the physical constants
!
      call setcon
!
!   Open the  .i  file
!
      call setsum(name,nci)
!
!   Open, check, load data from, and close, the  .csl  file
!
      call setcsla(name,ncore_not_used)
!
!   Get the remaining information
!
      call getsmd(name)
!
!   Get the eigenvectors
!
!      PRINT *, 'Block format?'
!      YES = GETYN ()
!      PRINT *

!      IF (YES) THEN
         call getmixblock(name,nci)
!      ELSE
!         IF (NCI.EQ.0) THEN
!            CALL GETMIXC(NAME)
!         ELSE
!            CALL GETMIXA(NAME)
!         ENDIF
!      ENDIF
!
!   Append a summary of the inputs to the  .sum  file
!
!      CALL STRSUM
!
!   Set up the table of logarithms of factorials
!
      call factt
!
!   Proceed with the RIS calculation
!
      call ris_cal(name)
!
!   Print completion message
!
      write (stderr_unit, *)
!      CALL STOPTIME (ncount1, 'RIS')
      write (stderr_unit, *) "RIS: Execution complete."
!
      stop
      end  program ris
