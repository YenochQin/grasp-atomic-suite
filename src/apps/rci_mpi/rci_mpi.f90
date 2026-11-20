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
      program rci_mpi
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
      use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_process_count, mpi_rank, &
         & processor_name, processor_name_length
      use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use calculation_state_files, only: read_state_name, orbital_input_path
      use yes_no_input, only: getyn
      use mpi_csf_block_input, only: open_csf_blocks_mpi
      use setdbg_i, only: setdbg
      use machine_constants, only: setmc
      use physical_constants, only: setcon
      use setsum_i, only: setsum
      use setres_i, only: setres
      use setmixmpi_i, only: setmixmpi
      use strsum_i, only: strsum
      use log_factorial_table, only: factt
      use matrix_i, only: matrix
      use genintbreit1wrap_i, only: genintbreit1wrap
      use genintbreit2wrap_i, only: genintbreit2wrap
      use genintrkwrap_i, only: genintrkwrap
      use mpi_utilities, only: startmpi2
      use mpi_utilities, only: stopmpi2
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nblk0 = 50
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, parameter :: nblk=50
      integer :: start_clock_count, finish_clock_count, ncount_rate, ncount_max
      integer :: nymduhmsm(8)
      integer :: core_subshell_count, ndum, max_orbital_two_j, nseconds
      logical :: yes
!      CHARACTER(LEN=128) :: NAME, TMPDIR, PERMDIR, ISOFILE
      character(len=128) :: startup_directory, state_name, tmpdir, permdir
      character(:), allocatable :: name, isotope_file_path, orbital_path
      character(len=8) :: idblk(nblk0)
      character :: chdate*8, chtime*10, chzone*5, str*8, msg*128
!-----------------------------------------------
!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================

      call startmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length, start_clock_count,       &
                           startup_directory, permdir, tmpdir, "RCI_MPI")


!=======================================================================
!  Get NDEF on node-0 and then send to all nodes
!=======================================================================

      if (mpi_rank == 0) then
         write (stderr_unit,*)
         write (stderr_unit, *) "Default settings? "
         yes = getyn()
         if (yes) then
            ndef = 0
!cjb fort.734
!cjb        write(734,'(A)') 'y            ! Default settings'
         else
            ndef = 1
!cjb fort.734
!cjb        write(734,'(A)') 'n            ! Default settings'
         end if
      end if
      call mpi_bcast (ndef,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
!     CALL MPI_Bcast (imethod,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


!=======================================================================
!  Get name of the state (used in files like <name>.c, <name>.s)
!  info used only on node-0
!=======================================================================

      name = ""
      isotope_file_path = ""
      orbital_path = ""
      if (mpi_rank == 0) then
         call read_state_name(state_name)
         open(unit=734,file=trim(state_name)//".clog",status="unknown")
         write(734,"(a)") "y            ! Default settings"
         write(734,"(a)") trim(state_name)
         isotope_file_path = trim(permdir)//"/isodata"
         name = trim(permdir)//"/"//trim(state_name)
         orbital_path = orbital_input_path(name)
      end if

!=======================================================================
   99 continue


      imcdf = 26        ! Unit for rci.res file
      iprerun = 0
!=======================================================================

!=======================================================================
!
!  Checks and settings... Mostly done in backyard.
!
!    SETDBG - Debug output control parameters - all set to false
!    SETMC - machine- and installation- dependent constants
!    SETCON - physical constants
!    SETSUM - open the summary file
!    cslhmpi - load header of the CSL file
!    SETRES - setup the .res  file, one for each node
!    SETMIXmpi - set up mixing coefficients file AND process user input
!    FACTT - table of logarithms of factorials setup
!    SETCSLmpi - open, check, load data from, and close the  .csl  file
!=======================================================================

      call setdbg
      call setmc
      call setcon
      if (mpi_rank == 0)  call setsum (name)
      call open_csf_blocks_mpi (name // "_descriptors.parquet", core_subshell_count, idblk)
      call setres (isotope_file_path, orbital_path, idblk)
      call setmixmpi (name, idblk)
      if (mpi_rank == 0)  call strsum
      call factt

!=======================================================================
!   Calculate all the needed Rk integrals - genintrkwrap calls genintrk
!   to compute the integrals (each processor allocates the same amount
!   of memory, but only computes part of the integrals); and then do the
!   "MPI_Allreduce" stuff.
!=======================================================================

      call genintrkwrap (mpi_rank, mpi_process_count, max_orbital_two_j)

!
!   If transverse interaction comput Breit integrals of type 1 and 2
!
      if (ltrans) then
         call genintbreit1wrap (mpi_rank, mpi_process_count, max_orbital_two_j)
         call genintbreit2wrap (mpi_rank, mpi_process_count, max_orbital_two_j)
      end if

!=======================================================================
!   Proceed with the CI calculation
!=======================================================================
      write (6, *) "Calling MATRIX..."

      call matrix (core_subshell_count, max_orbital_two_j)

      if (iprerun == 1) then
         iprerun = 2
         goto 99
      end if
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================

      call stopmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length,          &
                           start_clock_count, "RCI_MPI")

      stop
      end program rci_mpi
