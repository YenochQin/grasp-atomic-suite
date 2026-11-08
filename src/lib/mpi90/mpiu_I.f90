module mpiu_I
    implicit none
    interface
      subroutine startmpi (myid, nprocs, host, lenhost)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER, INTENT(OUT)  :: myid,  nprocs, lenhost
      CHARACTER(LEN=*), INTENT(OUT) :: host
      end subroutine startmpi
      subroutine mpix_startup (myid, nprocs, host, lenhost, ncount1, &
                               progname, description, infiledesc, &
                               outfiledesc, quiet_workers)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER :: myid, nprocs, lenhost, ncount1
      CHARACTER(LEN=*), INTENT(OUT) :: host
      CHARACTER(LEN=*), INTENT(IN) :: progname, description
      CHARACTER(LEN=*), INTENT(IN) :: infiledesc, outfiledesc
      LOGICAL, INTENT(IN) :: quiet_workers
      end subroutine mpix_startup
      subroutine startmpi2 (myid, nprocs, host, lenhost, ncount1, &
                           startdir, permdir, tmpdir, progname)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER, INTENT(OUT)  :: myid,  nprocs, lenhost
      INTEGER, INTENT(OUT)  :: ncount1
      CHARACTER(LEN=*), INTENT(OUT) :: host
      CHARACTER(LEN=*), INTENT(IN) :: progname
      CHARACTER(LEN=*), INTENT(OUT) :: startdir, permdir, tmpdir
      end subroutine startmpi2
      subroutine mpix_shutdown (myid, ncount1, progname)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER, INTENT(IN) :: myid, ncount1
      CHARACTER(LEN=*), INTENT(IN) :: progname
      end subroutine mpix_shutdown
      subroutine stopmpi (what, myid)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER, INTENT(IN)  ::  myid
      CHARACTER(LEN=*), INTENT(IN) :: WHAT
      end subroutine stopmpi
      subroutine stopmpi2 (myid, nprocs, host, lenhost, ncount1, &
                            progname)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER, INTENT(IN)  ::  myid,nprocs,lenhost
      integer              ::  ncount1
      character*(*), INTENT(IN) ::  host, progname
      end subroutine stopmpi2
      subroutine mpix_printmsg (msg, myid, nprocs)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER, INTENT(IN)  ::  myid, nprocs
      CHARACTER(LEN=*), INTENT(INOUT) :: msg
      end subroutine mpix_printmsg
      subroutine mpix_chkpt (myid, what)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER, INTENT(IN) :: myid
      CHARACTER(LEN=*), INTENT(IN) :: what
      end subroutine mpix_chkpt
      subroutine mpix_bytes (n, newType, ierr)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(out) :: newtype, ierr
      end subroutine mpix_bytes
      subroutine gdsummpi (x, n)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER, INTENT(IN)                         :: n
      REAL(DOUBLE), DIMENSION(1:n), INTENT(INOUT) ::   x
      end subroutine gdsummpi
      subroutine gdrsummpi_root (x, n)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER, INTENT(IN)                          :: n
      REAL(DOUBLE), DIMENSION(1:n), INTENT(INOUT) :: x
      end subroutine gdrsummpi_root
      subroutine gdmaxmpi_root (x, n)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER, INTENT(IN)                          :: n
      REAL(DOUBLE), DIMENSION(1:n), INTENT(INOUT) :: x
      end subroutine gdmaxmpi_root
      subroutine gisummpi (ix, n)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER, INTENT(IN)                         :: n
      INTEGER, DIMENSION(1:n), INTENT(INOUT) ::  ix
      end subroutine gisummpi
    end interface
end module
