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
      use default_c, only: double, inpci, ndef, ndump
      use blim_c, only: coeffcut1, coeffcut2, double, iprerun, ncsfpre
      use where_c
      use cons_c
      use def_c
      use hblock_c, only: csf_counts_by_block, double, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use blk_c
      use iccu_c, only: iccut, iccutblk
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      use terms_c, only: itab, jtab, nrows, ntab
      use mpi_c
      use decide_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i, only: getyn
      use cslhmpi_i, only: cslhmpi
      use setdbg_i, only: setdbg
      use setmc_i, only: setmc
      use setcon_i, only: setcon
      use setsum_i, only: setsum
      use setres_i, only: setres
      use setmixmpi_i, only: setmixmpi
      use strsum_i, only: strsum
      use factt_i, only: factt
      use matrix_i, only: matrix
      use genintbreit1wrap_i, only: genintbreit1wrap
      use genintbreit2wrap_i, only: genintbreit2wrap
      use genintrkwrap_i, only: genintrkwrap
      use mpiu_i, only: startmpi2
      use mpiu_i, only: stopmpi2
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nblk0 = 50
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, parameter :: nblk=50
      integer            :: lenperm, lentmp
      integer :: ncount1, ncount2, ncount_rate, ncount_max
      integer :: nymduhmsm(8)
      integer :: k, lenname, ncore, ndum, j2max, nseconds
      logical :: yes
      character (len = 3) :: idstring
!      CHARACTER(LEN=128) :: NAME, TMPDIR, PERMDIR, ISOFILE
      character(len=128) :: startdir, name, isofile, namesave,tmpdir, permdir
      character(len=8) :: idblk(nblk0)
      character :: chdate*8, chtime*10, chzone*5, str*8, msg*128
!-----------------------------------------------
!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================

      call startmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length, ncount1,       &
                           startdir, permdir, tmpdir, "RCI_MPI")
      write (idstring, "(I3.3)") mpi_rank

      lenperm = len_trim (permdir)
      lentmp = len_trim (tmpdir)

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

      lenname = 0     ! Otherwise it might cause trouble for myid!=0
      if (mpi_rank == 0) then
         do
            write (stderr_unit,*) "Name of state: "
            read (*,"(A)") name
            k = index (name," ")
            if (k > 1) exit
            write (stderr_unit,*) "Name may not start with a blank. redo..."
         end do
         namesave = name

         if (mpi_rank == 0) then
            open(unit=734,file=trim(name)//".clog",status="unknown")
            write(734,"(a)") "y            ! Default settings"
            write(734,"(a)") trim(name)
         end if

         !...Form the full name of the files used on node-0

         lenname = len_trim (name)
         isofile = permdir(1:lenperm) // "/isodata "
         name    = permdir(1:lenperm) // "/" // name(1:lenname) // " "
         lenname = lenperm + lenname + 1
!CFF     .. already open
!        open(unit=734,file=trim(name)//'.clog',status='unknown')
!        write(734,'(a)') 'y            ! Default settings'
!        write(734,'(a)') trim(namesave)
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
      call cslhmpi (name(1:lenname) // ".c", ncore, nblk0, idblk)
      call setres (isofile, name(1:lenname) // ".w", idblk)
      call setmixmpi (name, idblk)
      if (mpi_rank == 0)  call strsum
      call factt

!=======================================================================
!   Calculate all the needed Rk integrals - genintrkwrap calls genintrk
!   to compute the integrals (each processor allocates the same amount
!   of memory, but only computes part of the integrals); and then do the
!   "MPI_Allreduce" stuff.
!=======================================================================

      call genintrkwrap (mpi_rank, mpi_process_count, j2max)

!
!   If transverse interaction comput Breit integrals of type 1 and 2
!
      if (ltrans) then
         call genintbreit1wrap (mpi_rank, mpi_process_count, j2max)
         call genintbreit2wrap (mpi_rank, mpi_process_count, j2max)
      end if

!=======================================================================
!   Proceed with the CI calculation
!=======================================================================
      write (6, *) "Calling MATRIX..."

      call matrix (ncore, j2max)

      if (iprerun == 1) then
         iprerun = 2
         goto 99
      end if
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================

      call stopmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length,          &
                           ncount1, "RCI_MPI")

      stop
      end program rci_mpi
