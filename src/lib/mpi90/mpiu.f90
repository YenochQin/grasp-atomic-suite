! Collection of FORTRAN subroutines for general purpose MPI operations.
! Mostly IMPLICIT NONE
! MPI calls use the implementation's explicit module interfaces.
!***********************************************************************
      subroutine startmpi (myid, nprocs, host, lenhost)
! This subroutine starts MPI to get the process id (myid), the number
! of processes (nprocs) and the status of enrollment (ierr)

!***********************************************************************
!  Modified by Charlotte F. Fischer   10/10/2017
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      use mpi
      implicit none
      integer, intent(out)  :: myid,  nprocs, lenhost
      character(len=*), intent(out) :: host
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      integer ierr

      call mpi_init (ierr)
      call mpi_comm_rank (mpi_comm_world, myid, ierr)
      call mpi_comm_size (mpi_comm_world, nprocs, ierr)
      call mpi_get_processor_name (host, lenhost, ierr)
      if (ierr .ne. mpi_success) print *, ' ierr=', ierr, ' id=', myid

      return
      end

!***********************************************************************
      subroutine mpix_startup (myid, nprocs, host, lenhost, ncount1, &
                               progname, description, infiledesc, &
                               outfiledesc, quiet_workers)
! Shared-filesystem-safe MPI startup for the post-processing programs.
! Unlike startmpi2, this variant does not change the working directory.
!***********************************************************************
      use mpi
      use mpiu_i, only: startmpi
      implicit none
      integer :: myid, nprocs, lenhost, ncount1
      character(len=*), intent(out) :: host
      character(len=*), intent(in) :: progname, description
      character(len=*), intent(in) :: infiledesc, outfiledesc
      logical, intent(in) :: quiet_workers
      integer :: ncount_rate, ncount_max, iouterr

      call startmpi (myid, nprocs, host, lenhost)
      call system_clock (ncount1, ncount_rate, ncount_max)

      if (myid .eq. 0) then
         write (6, *)
         write (6, *) trim(progname)
         write (6, *) trim(description)
         write (6, *) 'Input files:  ' // trim(infiledesc)
         write (6, *) 'Output files: ' // trim(outfiledesc)
         write (6, *) 'Running on ', nprocs, ' MPI ranks'
      else if (quiet_workers) then
         open (unit=6, file='/dev/null', status='OLD', action='WRITE', &
              iostat=iouterr)
      endif

      return
      end

!***********************************************************************
      subroutine startmpi2 (myid, nprocs, host, lenhost, ncount1, &
                           startdir, permdir, tmpdir, progname)
! Calls startmpi to get mpi environment;
! Calls cpath to get various paths;
! Calls DATE_AND_TIME to get date, time, zone;
! Calls mpi_printmsg to print some of these info to screen.

!***********************************************************************
!  Modified by Charlotte F. Fischer   10/10/2017

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      use cpath_i, only: cpath
      use mpiu_i, only: mpix_printmsg
      use mpiu_i, only: startmpi
      implicit none
      integer, intent(out)  :: myid,  nprocs, lenhost
      integer, intent(out)  :: ncount1
      character(len=*), intent(out) :: host
      character(len=*), intent(in) :: progname
      character(len=*), intent(out) :: startdir, permdir, tmpdir
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      character idstring*3
      character chdate*8, chtime*10, chzone*5
               !ccyymmdd  hhmmss.sss  Shhmm
      integer  nymduhmsm(8)
               !Year Month Day Universal Hour Minute Sesond Millisecond
      integer lenstart, lenperm, lentmp, ncount_rate, ncount_max
! For printing
      character(len=80) :: msg

!=======================================================================
!  Get processor info: myid, nprocs, host name; and print
!=======================================================================

      call startmpi (myid, nprocs, host, lenhost)
      write (idstring, '(I3.3)') myid
      if (myid .eq. 0) then
         print *, '===================================================='
         print *, '       ', progname, ': Execution Begins ...'
         print *, '===================================================='
         print *,        'Participating nodes:'
       endif
      msg = '  Host: ' // host(1:lenhost) // '    ID: ' // idstring
      call mpix_printmsg (msg, myid, nprocs)

!=======================================================================
!  Get date, time, zone and print
!=======================================================================

      call date_and_time (chdate, chtime, chzone, nymduhmsm)
      if (myid .eq. 0) then
         print *
         print *, 'Date and Time:'
      endif
      msg = '  ' // host(1:lenhost) // ': ' // &
            '  Date: ' // chdate // &
            '  Time: ' // chtime // &
            '  Zone: ' // chzone
      call mpix_printmsg (msg, myid, nprocs)

!=======================================================================
!  Set up local working dir and go there
!     tmpdir  - local working dir of the node. mcpXXX files are there
!     permdir - for I/O specific to node-0.
!=======================================================================

      call cpath (startdir, permdir, tmpdir)

      lenstart = len_trim (startdir)
      lenperm = len_trim (permdir)
      lentmp = len_trim (tmpdir)

      if (myid .eq. 0) then
         print *
         print *, 'Start Dir:'
      endif
      msg = '  ' // host(1:lenhost) // ': ' // startdir(1:lenstart)
      call mpix_printmsg (msg, myid, nprocs)

      if (myid .eq. 0) then
         print *
         print *, 'Serial I/O Dir (node-0 only):'
         print *, '  ' // host(1:lenhost) // ': ' // permdir(1:lenperm)
      endif

      if (myid .eq. 0) then
         print *
         print *, 'Work Dir (Parallel I/O):'
      endif
      msg = '  ' // host(1:lenhost) // ': ' // tmpdir(1:lentmp)
      call mpix_printmsg (msg, myid, nprocs)

!=======================================================================
!  Start timing - Record the wall clock
!=======================================================================

      call system_clock (ncount1, ncount_rate, ncount_max)
      return
      end

!***********************************************************************
      subroutine mpix_shutdown (myid, ncount1, progname)
! Shared-filesystem-safe MPI shutdown companion to mpix_startup.
!***********************************************************************
      use mpi
      implicit none
      integer, intent(in) :: myid, ncount1
      character(len=*), intent(in) :: progname
      integer :: ncount2, ncount_rate, ncount_max, nseconds, ierr
      character(len=8) :: chdate
      character(len=10) :: chtime
      character(len=5) :: chzone
      integer, dimension(8) :: nymduhmsm

      call system_clock (ncount2, ncount_rate, ncount_max)
      ncount2 = ncount2 - ncount1
      nseconds = ncount2 / ncount_rate
      call date_and_time (chdate, chtime, chzone, nymduhmsm)

      if (myid .eq. 0) then
         write (6, *)
         write (6, *) trim(progname) // ': Execution complete.'
         write (6, *) 'Wall time: ', nseconds, ' seconds'
         write (6, *) 'Finish date: ', chdate, ' time: ', chtime, &
                      ' zone: ', chzone
      endif

      call mpi_finalize (ierr)

      return
      end

!***********************************************************************
      subroutine stopmpi (what, myid)

! This subroutine stops MPI, but before doing so it issues an error
! message to node-0 about what went wrong and where it happened.
!   what - string, type of error (most likely a subroutine name)
!   myid - The id of the PE where error happened

!***********************************************************************
!  Modified by Charlotte F. Fischer   10/10/2017

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      use mpi
      implicit none
      integer, intent(in)  ::  myid
      character(len=*), intent(in) :: what
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer ierr

      print *, 'mpi stopped by node-', myid, ' from ', what
      call mpi_barrier (mpi_comm_world,ierr)
      call mpi_finalize (ierr)

      stop
      end

!***********************************************************************
      subroutine stopmpi2 (myid, nprocs, host, lenhost, ncount1, &
                            progname)

! Calls DATE_AND_TIME to get date, time, zone;
! Calls mpi_printmsg to print some of these info to screen.

!***********************************************************************
!  Modified by Charlotte F. Fischer   10/10/2017

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      use mpi
      use mpiu_i, only: mpix_printmsg
      use mpiu_i, only: stopmpi
      implicit none
      integer, intent(in)  ::  myid,nprocs,lenhost
      integer              ::  ncount1
      character*(*), intent(in) ::  host, progname
! Things for timing
      integer   ncount2, ncount_rate, ncount_max, nseconds
      character chdate*8, chtime*10, chzone*5
               !ccyymmdd  hhmmss.sss  Shhmm
      integer  nymduhmsm(8)
               !Year Month Day Universal Hour Minute Sesond Millisecond

! For printing
      character str2nds*8, msg*80

      integer ierr

!=======================================================================
!  Get processor info: myid, nprocs, host name; and print
!=======================================================================

      if (myid .eq. 0) then
         print *, '===================================================='
         print *, '       ', progname, ': Execution Finished ...'
         print *, '===================================================='
         print *,        'Wall time:'
       endif

      call system_clock (ncount2, ncount_rate, ncount_max)
      ncount2 = ncount2 - ncount1
      nseconds = ncount2 / ncount_rate
      write (str2nds, '(i8)') nseconds
      msg = str2nds // ' seconds on ' // host(1:lenhost)
      call mpix_printmsg (msg, myid, nprocs)

      if (myid .eq. 0) then
         print *
         print *, 'Finish Date and Time:'
      endif

      call date_and_time (chdate, chtime, chzone, nymduhmsm)

      msg = '  ' // host(1:lenhost) // ': ' // &
            '  Date: ' // chdate // &
            '  Time: ' // chtime // &
            '  Zone: ' // chzone
      call mpix_printmsg (msg, myid, nprocs)

      if (myid .eq. 0) print *

      call mpi_barrier (mpi_comm_world,ierr)
      call stopmpi (progname // ': Execution complete.', myid)
      return
      end

!***********************************************************************
      subroutine mpix_printmsg (msg, myid, nprocs)

! Displays on node-0's screen info from all nodes including node 0.

!***********************************************************************
!  Modified by Charlotte F. Fischer   10/10/2017

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      use mpi
      implicit none
      integer, intent(in)  ::  myid, nprocs
      character(len=*), intent(inout) :: msg

!-----------------------------------------------
!   L o c a l  V a r i a b l e s
!-----------------------------------------------
      integer :: inid, istat(mpi_status_size), ierr, msglength

      msglength = len_trim (msg)

      if (myid .ne. 0) then
         call mpi_send (msglength, 1, mpi_integer, 0, myid, &
                        mpi_comm_world, ierr)   ! Send nsgLength
         call mpi_send (msg, msglength, mpi_character, 0, myid+nprocs, &
                        mpi_comm_world, ierr)   ! Send msg
      else
         print *, msg(1:msglength)      ! msg from node 0 itself
         do inid = 1, nprocs - 1
            call mpi_recv (msglength, 1, mpi_integer, inid, &
                           inid, mpi_comm_world, istat, ierr)
            call mpi_recv (msg, msglength, mpi_character, inid, &
                           inid+nprocs, mpi_comm_world, istat, ierr)
            print *, msg(1:msglength)
         enddo
      endif

      return
      end

!***********************************************************************
      subroutine mpix_chkpt (myid, what)
!     To set a check-point in mpi program
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      use mpi
      implicit none
      integer, intent(in) :: myid
      character(len=*), intent(in) :: what
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ierr
      character(len=1) :: chtmp
!-----------------------------------------------
!
      if (myid .eq. 0) then
         print *, what, 'hit any key to continue'
         read (*, '(a)') chtmp
      endif

      call mpi_barrier (mpi_comm_world, ierr)

      return
      end

!***********************************************************************
      subroutine mpix_bytes (n, newtype, ierr)

!  Constructs new mpi data type newType of n-bytes long
!***********************************************************************
!  Modified by Charlotte F. Fischer  10/10/2017
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      use mpi
      implicit none
      integer, intent(in) :: n
      integer, intent(out) :: newtype, ierr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer ier0

      call mpi_type_contiguous (n, mpi_byte, newtype, ier0)
      call mpi_type_commit (newtype, ierr)

      ierr = ierr + ier0

      return
      end

!***********************************************************************
      subroutine gdsummpi (x, n)

!     Sum x over all nodes and in y, then copy y to  x
!***********************************************************************
      use vast_kind_param, only:  double
      use mpi
      use suite_linear_algebra, only: dcopy
      use dinit_i, only: dinit
      implicit none

      integer, intent(in)                         :: n
      real(double), dimension(1:n), intent(inout) ::   x

      integer                                     :: ierr
      real(double), dimension(1:n)                :: y

      call dinit (n, 0.d0, y, 1)
      call mpi_allreduce (x, y, n, mpi_double_precision, &
                                   mpi_sum, mpi_comm_world, ierr)
      call dcopy (n, y, 1, x, 1)   ! copy y to x

      return
      end

!***********************************************************************
      subroutine gdrsummpi_root (x, n)

!     Sum x onto node-0 only; non-root ranks do not receive the result.
!***********************************************************************
      use vast_kind_param, only:  double
      use mpi_c, only: myid
      use mpi
      use suite_linear_algebra, only: dcopy
      use dinit_i, only: dinit
      implicit none

      integer, intent(in)                          :: n
      real(double), dimension(1:n), intent(inout) :: x

      integer                                      :: ierr
      real(double), dimension(1:n)                 :: y

      if (myid .eq. 0) then
         call dinit (n, 0.d0, y, 1)
         call mpi_reduce (x, y, n, mpi_double_precision, mpi_sum, 0, &
                          mpi_comm_world, ierr)
         call dcopy (n, y, 1, x, 1)
      else
         call mpi_reduce (x, y, n, mpi_double_precision, mpi_sum, 0, &
                          mpi_comm_world, ierr)
      endif

      return
      end

!***********************************************************************
      subroutine gdmaxmpi_root (x, n)

!     Max x onto node-0 only; non-root ranks do not receive the result.
!***********************************************************************
      use vast_kind_param, only:  double
      use mpi_c, only: myid
      use mpi
      use suite_linear_algebra, only: dcopy
      implicit none

      integer, intent(in)                          :: n
      real(double), dimension(1:n), intent(inout) :: x

      integer                                      :: ierr
      real(double), dimension(1:n)                 :: y

      if (myid .eq. 0) then
         call mpi_reduce (x, y, n, mpi_double_precision, mpi_max, 0, &
                          mpi_comm_world, ierr)
         call dcopy (n, y, 1, x, 1)
      else
         call mpi_reduce (x, y, n, mpi_double_precision, mpi_max, 0, &
                          mpi_comm_world, ierr)
      endif

      return
      end

!***********************************************************************
      subroutine gisummpi (ix, n)
!     Sum x over all nodes and in iy, then copy iy to  ix
!***********************************************************************
      use mpi
      use icopy_i, only: icopy
      implicit none

      integer, intent(in)                         :: n
      integer, dimension(1:n), intent(inout) ::  ix

      integer                                     :: ierr, i
      integer, dimension(1:n)                ::  iy

      iy = 0

      call mpi_allreduce (ix, iy, n, mpi_integer, mpi_sum, &
                                     mpi_comm_world, ierr)
      call icopy (n, iy, 1, ix, 1)

      return
      end
