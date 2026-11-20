

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
      use calculation_input_defaults, only: inpci, ndef, ndump
      use numeric_constants, only: eleven, eps, four, half, one, seven, ten, tenth, three, two, zero
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use rwfnestimate_input, only: radial_input_options, parse_radial_options, load_radial_header
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use calculation_state_files, only: check_orbital_output
      use calculation_state_files, only: state_files, state_name_length, select_state, files_for_state, orbital_input_path
      use setdbg_i, only: setdbg
      use machine_constants, only: setmc
      use physical_constants, only: setcon
      use setsum_i, only: setsum
      use screenpar_i, only: screenpar
      use getinf_i, only: getinf
      use strsum_i, only: strsum
      use genrwf_i, only: genrwf
      use radial_orbital_orthogonalization, only: orthsc
      use wrtrwf_i, only: wrtrwf
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      type(state_files) :: files
      character(state_name_length) :: state
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
      write (stderr_unit, *) "Input files: isodata, descriptor header TOML, optional rwfn file"
      call select_state(state, "Default settings ?", yes)
      files = files_for_state(state)
      call check_orbital_output(files%initial_orbitals, writing_initial=.true.)
      if (input_options%filename == "") input_options%filename = files%header
      block
         character(:), allocatable :: basename
         basename = input_options%filename(scan(input_options%filename, "/"//achar(92), back=.true.)+1:)
         if (basename == files%initial_orbitals .or. basename == files%estimate_debug .or. &
             basename == files%estimate_summary) error stop "Header uses a reserved output name"
      end block
      write (stderr_unit, *) "Selected input: ", input_options%filename
      write (stderr_unit, *) "Output file: ", files%initial_orbitals


!
      if (yes) then
         ndef = 0
      else
         ndef = 1
      end if

!
!   Determine if there is to be any debug printout; this will be
!   made on the  .dbg  file
!
      call setdbg(files%estimate_debug)
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
      if (ndef /= 0) call setsum(files%estimate_summary)
!
!   Load orbital metadata from the selected rCSFs header.
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
      if (state /= "") files%orbitals = orbital_input_path(state)
      call genrwf(files%orbitals)
!
!   Orthogonalize the radial orbitals
!
      call orthsc
!
!   Write the subshell radial wavefunctions out
!
      call wrtrwf(files%initial_orbitals)
!
!   Print completion message
!
      write (stderr_unit, *) "RWFNESTIMATE: Execution complete."
!
      stop
      end program erwf
