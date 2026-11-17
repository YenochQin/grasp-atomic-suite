module mpi_runtime_i
    implicit none (type, external)
    private
    public :: startmpi
    public :: mpix_bytes
    public :: stopmpi
contains
      subroutine startmpi (myid, nprocs, host, lenhost)
! This subroutine starts MPI to get the process id (myid), the number
! of processes (nprocs) and the status of enrollment (ierr)

!***********************************************************************
!  Modified by Charlotte F. Fischer   10/10/2017
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      use mpi, only: mpi_comm_rank, mpi_comm_size, mpi_comm_world, mpi_get_processor_name, mpi_init, mpi_success
      implicit none (type, external)
      integer, intent(out)  :: myid,  nprocs, lenhost
      character(len=*), intent(out) :: host
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      integer :: ierr

      call mpi_init (ierr)
      call mpi_comm_rank (mpi_comm_world, myid, ierr)
      call mpi_comm_size (mpi_comm_world, nprocs, ierr)
      call mpi_get_processor_name (host, lenhost, ierr)
      if (ierr /= mpi_success) print *, " ierr=", ierr, " id=", myid

      return
      end subroutine startmpi
      subroutine mpix_bytes (n, newtype, ierr)

!  Constructs new mpi data type newType of n-bytes long
!***********************************************************************
!  Modified by Charlotte F. Fischer  10/10/2017
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      use mpi, only: mpi_byte, mpi_type_commit, mpi_type_contiguous
      implicit none (type, external)
      integer, intent(in) :: n
      integer, intent(out) :: newtype, ierr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ier0

      call mpi_type_contiguous (n, mpi_byte, newtype, ier0)
      call mpi_type_commit (newtype, ierr)

      ierr = ierr + ier0

      return
      end subroutine mpix_bytes
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
      use mpi, only: mpi_barrier, mpi_comm_world, mpi_finalize
      implicit none (type, external)
      integer, intent(in)  ::  myid
      character(len=*), intent(in) :: what
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ierr

      print *, "mpi stopped by node-", myid, " from ", what
      call mpi_barrier (mpi_comm_world,ierr)
      call mpi_finalize (ierr)

      stop
      end subroutine stopmpi
end module mpi_runtime_i
