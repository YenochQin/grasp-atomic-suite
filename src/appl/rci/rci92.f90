!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**             ******    *****   ****   *****    *****              ***
!**             **   **  **   **   **   **   **  **   **             ***
!**             **   **  **        **   **   **       **             ***
!**             ******   **        **    *****       **              ***
!**             **  **   **        **      **       **               ***
!**             **   **  **   **   **     **      **                 ***
!**             **   **   *****   ****   **      *******             ***
!**                                                                  ***
!**          Relativistic Configuration-Interaction Program          ***
!**                                                                  ***
!**   This program is a derivative of GRASP2 (F. A. Parpia, I. P.    ***
!**   Grant, and C. F. Fischer, 1990).                               ***
!**                                                                  ***
!**                            GRASP92                               ***
!**          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
!**                                                                  ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
      program rci92
!                                                                      *
!   Entry routine for RCI92. Controls the entire computation.          *
!                                                                      *
!   Call(s) to: [LIB92]: SETMC, SETCON.                                *
!               [RCI92]: CHKPLT, MATRIX, SETCSL, SETDBG, SETMIX,       *
!                        SETRES, SETSUM, STRSUM.                       *
!               [NJGRAF]: FACTT.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Oct 1992   *
!   Updated by Xinghong He                Last revision: 23 Jun 1998   *
!   Modified by Gediminas Gaigalas for new spin-angular integration.   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use calculation_input_defaults, only: inpci, ndef, ndump
      use where_c, only: imcdf
      use rci_prerun_settings, only: coeffcut1, coeffcut2, iprerun, ncsfpre
      use numeric_constants, only: eleven, eps, four, half, one, seven, ten, tenth, three, two, zero
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use transition_basis_blocks, only: ncfblk, nevinblk, ncfinblk, two_j, nblock, idxblk, nelectot, ncftot, nwtot, nvectot, &
          nvecsiztot, nblock1, nblocki, nblockf, ncfi, ncff
      use csf_block_cutoffs, only: iccut, iccutblk
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use subshell_term_table, only: itab, jtab, nrows, ntab
      use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use calculation_state_files, only: read_state_name, orbital_input_path
      use yes_no_input, only: getyn
      use setdbg_i, only: setdbg
      use machine_constants, only: setmc
      use physical_constants, only: setcon
      use setsum_i, only: setsum
      use setcsl_i, only: setcsl
      use setres_i, only: setres
      use setmix_i, only: setmix
      use strsum_i, only: strsum
      use log_factorial_table, only: factt
      use matrix_i, only: matrix
      use genintbreit1_i, only: genintbreit1
      use genintbreit2_i, only: genintbreit2
      use genintrk_i, only: genintrk
      use run_timer_start, only: starttime
      use run_timer_finish, only: stoptime
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nblk0 = 50
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, parameter :: nblk=50
      integer :: ncount1, ncount2, ncount_rate, ncount_max
      integer :: nymduhmsm(8)
      integer :: myid, nprocs, lenname, ncore, ndum, j2max, nseconds
      logical :: yes
      character(len=128) :: name, tmpdir, permdir, isofile
      character(len=8) :: idblk(nblk0)
      character :: chdate*8, chtime*10, chzone*5, str*8, msg*128
!-----------------------------------------------
!

      imcdf = 26                                 ! Unit for rci.res file
      iprerun = 0
      myid = 0
      nprocs = 1

      write(*,*)
      write(*,*) "RCI"
      write(*,*) "This is the configuration interaction program "
      write(*,*) "Input file:  isodata, name_descriptors.parquet, name_header.toml, name.wfn"
      write(*,*) "Outputfiles: name.cimix, name.csum, name.clog "
      write(*,*) "             rci.res (can be used for restart)"
      write(*,*)


!
! Start timing
!
      call starttime (ncount1, "RCI")
!
! Get NDEF
!
!      WRITE (ISTDE, *) 'RCI2: Execution begins ...'
!      WRITE (ISTDE, *)
      write (stderr_unit, "(A)") "Default settings? "
      yes = getyn()
      if (yes) then
         ndef = 0
      else
         ndef = 1
      end if
!
! Get name of the state (used in files like <name>.c, <name>.s)
!
      call read_state_name(name)

! Now the name of the state is known, open the log file

      open(unit=734, file=trim(name)//".clog",status="unknown")
      write(734,"(a)") "y            ! Default settings"
      write(734,"(a)") trim(name)

!         ...Form the full name of the files used on node-0

      lenname = len_trim(name)
      isofile = "isodata"
!      WRITE (6, *) 'isofile = ', ISOFILE(1:LEN_TRIM(ISOFILE))
!      WRITE (6, *) 'name = ', NAME(1:LEN_TRIM(NAME))

   99 continue

!
! In SETDBG of this version all control logicals are set to
! false thus no debug output will be made
!
!      WRITE (6, *) 'Calling SETDBG...'
      call setdbg
!
! Perform machine- and installation-dependent setup
!
!      WRITE (6, *) 'Calling SETMC...'
      call setmc
!
! Set up the physical constants
!
!      WRITE (6, *) 'Calling SETCON...'
      call setcon
!
! Open summary file
!
!      WRITE (6, *) 'Calling SETSUM...'
      call setsum (name)

!      WRITE (6, *) 'Calling setcsl...'
      call setcsl (name(1:lenname)//"_descriptors.parquet", ncore, nblk0, idblk)
!
! Set up the  .res  file; determine if this is a restart.
!
!      WRITE (6, *) 'Calling SETRES...'
      call setres (isofile, orbital_input_path(name), idblk)
!
!   Open the  .mix  file; determine the eigenpairs required
!
!      WRITE (6, *) 'Calling SETMIX...'
      call setmix (name, idblk)
!
!   Append a summary of the inputs to the  .sum  file
!
      write (6, *) "Calling STRSUM..."
      call strsum
!
!   Set up the table of logarithms of factorials
!
      write (6, *) "Calling FACTT..."
      call factt
!
!   Calculate all the needed Rk integrals
!
      write (6, *) "Calling GENINTRK..."
      call genintrk (myid, nprocs, ndum, j2max)
!
!   If transverse interaction comput Breit integrals of type 1 and 2
!
      if (ltrans) then
         print *, "Calling GENINTBREIT1..."
         call genintbreit1 (myid, nprocs, ndum, j2max)
         print *, "Calling GENINTBREIT2..."
         call genintbreit2 (myid, nprocs, ndum, j2max)
      end if
!
!   Proceed with the CI calculation
!
      write (6, *) "Calling MATRIX..."

      call matrix (ncore, j2max)

      if (iprerun == 1) then
         iprerun = 2
         goto 99
      end if

      if (myid == 0) then
         write (6, *)
         write (6, *)
         write (6, *) "Finish time, Statistics"
         write (6, *)
      end if
      close(734)
      call stoptime (ncount1, "RCI")
!      CALL SYSTEM_CLOCK (NCOUNT2, NCOUNT_RATE, NCOUNT_MAX)
!      NCOUNT2 = NCOUNT2 - NCOUNT1
!      NSECONDS = NCOUNT2/NCOUNT_RATE
!      WRITE (STR, '(I8)') NSECONDS
!      MSG = STR//' seconds '
!      WRITE (6, *) MSG
!
!      CALL DATE_AND_TIME (CHDATE, CHTIME, CHZONE, NYMDUHMSM)
!
!      MSG = ' Date: '//CHDATE//' Time: '//CHTIME//' Zone: '//CHZONE
!      WRITE (6, *) MSG
!
!   Print completion message
!
!      WRITE (6, *) 'RCI2: Execution complete.'
!
      stop
      end program rci92
