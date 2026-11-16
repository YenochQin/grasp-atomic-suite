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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use default_c
      use blim_c
      use where_c
      use cons_c
      use def_c
      use hblock_c
      use blk_c
      use iccu_c
      use iounit_c
      use terms_c
      use mpi_c
      use decide_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use cslhmpi_i
      use setdbg_i
      use setmc_i
      use setcon_i
      use setsum_i
      use setres_i
      use setmixmpi_i
      use strsum_i
      use factt_i
      use matrix_i
      use genintbreit1wrap_i, only: genintbreit1wrap
      use genintbreit2wrap_i, only: genintbreit2wrap
      use genintrkwrap_i, only: genintrkwrap
      use mpiu_i, only: startmpi2
      use mpiu_i, only: stopmpi2
      implicit none
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
      integer, dimension(8) :: nymduhmsm
      integer :: k, lenname, ncore, ndum, j2max, nseconds
      logical :: yes
      character (len = 3) :: idstring
!      CHARACTER(LEN=128) :: NAME, TMPDIR, PERMDIR, ISOFILE
      character(len=128) :: startdir, name, isofile, namesave,tmpdir, permdir
      character(len=8), dimension(nblk0) :: idblk
      character :: chdate*8, chtime*10, chzone*5, str*8, msg*128
!-----------------------------------------------
!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================

      call startmpi2 (myid, nprocs, host, lenhost, ncount1,       &
                           startdir, permdir, tmpdir, 'RCI_MPI')
      write (idstring, '(I3.3)') myid

      lenperm = len_trim (permdir)
      lentmp = len_trim (tmpdir)

!=======================================================================
!  Get NDEF on node-0 and then send to all nodes
!=======================================================================

      if (myid .eq. 0) then
         write (istde,*)
         write (istde, *) 'Default settings? '
         yes = getyn()
         if (yes) then
            ndef = 0
!cjb fort.734
!cjb        write(734,'(A)') 'y            ! Default settings'
         else
            ndef = 1
!cjb fort.734
!cjb        write(734,'(A)') 'n            ! Default settings'
         endif
      endif
      call mpi_bcast (ndef,1,mpi_integer,0,mpi_comm_world,ierr)
!     CALL MPI_Bcast (imethod,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


!=======================================================================
!  Get name of the state (used in files like <name>.c, <name>.s)
!  info used only on node-0
!=======================================================================

      lenname = 0     ! Otherwise it might cause trouble for myid!=0
      if (myid .eq. 0) then
         do
            write (istde,*) 'Name of state: '
            read (*,'(A)') name
            k = index (name,' ')
            if (k .gt. 1) exit
            write (istde,*) 'Name may not start with a blank. redo...'
         enddo
         namesave = name

         if (myid .eq. 0) then
            open(unit=734,file=trim(name)//'.clog',status='unknown')
            write(734,'(a)') 'y            ! Default settings'
            write(734,'(a)') trim(name)
         endif

         !...Form the full name of the files used on node-0

         lenname = len_trim (name)
         isofile = permdir(1:lenperm) // '/isodata '
         name    = permdir(1:lenperm) // '/' // name(1:lenname) // ' '
         lenname = lenperm + lenname + 1
!CFF     .. already open
!        open(unit=734,file=trim(name)//'.clog',status='unknown')
!        write(734,'(a)') 'y            ! Default settings'
!        write(734,'(a)') trim(namesave)
      endif

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
      if (myid .eq. 0)  call setsum (name)
      call cslhmpi (name(1:lenname) // '.c', ncore, nblk0, idblk)
      call setres (isofile, name(1:lenname) // '.w', idblk)
      call setmixmpi (name, idblk)
      if (myid .eq. 0)  call strsum
      call factt

!=======================================================================
!   Calculate all the needed Rk integrals - genintrkwrap calls genintrk
!   to compute the integrals (each processor allocates the same amount
!   of memory, but only computes part of the integrals); and then do the
!   "MPI_Allreduce" stuff.
!=======================================================================

      call genintrkwrap (myid, nprocs, j2max)

!
!   If transverse interaction comput Breit integrals of type 1 and 2
!
      if (ltrans) then
         call genintbreit1wrap (myid, nprocs, j2max)
         call genintbreit2wrap (myid, nprocs, j2max)
      end if

!=======================================================================
!   Proceed with the CI calculation
!=======================================================================
      write (6, *) 'Calling MATRIX...'

      call matrix (ncore, j2max)

      if (iprerun == 1) then
         iprerun = 2
         go to 99
      endif
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================

      call stopmpi2 (myid, nprocs, host, lenhost,          &
                           ncount1, 'RCI_MPI')

      stop
      end program rci_mpi
