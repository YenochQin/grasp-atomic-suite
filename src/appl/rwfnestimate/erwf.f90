

!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**               *******  *******   **    **  *******               ***
!**               **       **    **  **    **  **                    ***
!**               **       **    **  **    **  **                    ***
!**               *****    *******   ** ** **  *****                 ***
!**               **       **  **    ** ** **  **                    ***
!**               **       **   **   ** ** **  **                    ***
!**               *******  **    **   **  **   **                    ***
!**                                                                  ***
!**   Program to create a  .rwf  file for RSCF92 by merging one or   ***
!**   more existing  .rwf  files or generating Thomas-Fermi or hy-   ***
!**   drogenic functions.                                            ***
!**                                                                  ***
!**                              GRASP92                             ***
!**           F. A. Parpia, C. F. Fischer, and I. P. Grant           ***
!**                                                                  ***
!**                                                                  ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
      program erwf
!                                                                      *
!   Entry routine for RCI92. Controls the entire computation.          *
!                                                                      *
!   Call(s) to: [LIB92]: SETMC, SETCON.                                *
!               [RCI92]: CHKPLT, MATRIX, SETCSL, SETDBG, SETMIX,       *
!                        SETRES, SETSUM, STRSUM.                       *
!               [NJGRAF]: FACTT.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use default_c, only: double, inpci, ndef, ndump
      use cons_c, only: double, eleven, eps, four, half, one, seven, ten, tenth, three, two, zero
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      use rwfnestimate_input, only: radial_input_options, parse_radial_options, load_radial_header
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i, only: getyn
      use setdbg_i, only: setdbg
      use setmc_i, only: setmc
      use setcon_i, only: setcon
      use setsum_i, only: setsum
      use screenpar_i, only: screenpar
      use getinf_i, only: getinf
      use strsum_i, only: strsum
      use genrwf_i, only: genrwf
      use orthsc_i, only: orthsc
      use wrtrwf_i, only: wrtrwf
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncore
      logical :: yes
      type(radial_input_options) :: input_options
!-----------------------------------------------

!   Startup message
!
      call parse_radial_options(input_options)
      write (stderr_unit, *) "RWFNESTIMATE"
      write (stderr_unit, *) "This program estimates radial wave functions"
      write (stderr_unit, *) "for orbitals"
      write (stderr_unit, *) "Input files: isodata, selected CSF text/header, optional rwfn file"
      write (stderr_unit, *) "Selected input: ", input_options%filename
      write (stderr_unit, *) "Output file: rwfn.inp"


!
      write (stderr_unit, *) "Default settings ?"
      yes = getyn()
      if (yes) then
         ndef = 0
      else
         ndef = 1
      end if

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
!   Open the  .sum  file
!
      if (ndef /= 0) call setsum
!
!   Load orbital metadata from the selected CSF text or rCSFs header.
!
      call load_radial_header(input_options, ncore)
!
!   Hydrogenic screen parameters for all orbitals
!
      call screenpar (ncore)
!
!   Determine other relevant information
!
      call getinf
!
!   Write the first part of the  .sum  file
!
      if (ndef /= 0) call strsum
!
!   Generate the subshell radial wavefunctions
!
      call genrwf
!
!   Orthogonalize the radial orbitals
!
      call orthsc
!
!   Write the subshell radial wavefunctions out
!
      call wrtrwf
!
!   Print completion message
!
      write (stderr_unit, *) "RWFNESTIMATE: Execution complete."
!
      stop
      end program erwf
