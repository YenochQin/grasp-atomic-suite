      MODULE mpi_C
      USE vast_kind_param, ONLY:  DOUBLE
!...Created by Pacific-Sierra Research 77to90  4.3E  06:23:52  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
! Use the MPI implementation's explicit Fortran interfaces. Integer handles
! remain the public convention of this legacy module; choice buffers are
! described by the MPI module rather than unchecked external declarations.
      USE mpi
      IMPLICIT NONE

      character*(MPI_MAX_PROCESSOR_NAME) host
      INTEGER myid, nprocs, lenhost
      INTEGER ierr, ierrtotal
      INTEGER istat(MPI_STATUS_SIZE)

      INTERFACE
         SUBROUTINE startmpi(myid, nprocs, host, lenhost)
            INTEGER, INTENT(OUT) :: myid, nprocs, lenhost
            CHARACTER(*), INTENT(OUT) :: host
         END SUBROUTINE startmpi
         SUBROUTINE mpix_bytes(n, newtype, ierr)
            INTEGER, INTENT(IN) :: n
            INTEGER, INTENT(OUT) :: newtype, ierr
         END SUBROUTINE mpix_bytes
      END INTERFACE

      END MODULE mpi_C
