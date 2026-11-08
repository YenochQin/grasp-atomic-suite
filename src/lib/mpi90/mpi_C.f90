      MODULE mpi_C
      USE vast_kind_param, ONLY:  DOUBLE
!...Created by Pacific-Sierra Research 77to90  4.3E  06:23:52  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
! Use the MPI implementation's explicit Fortran interfaces. Integer handles
! remain the public convention of this legacy module; choice buffers are
! described by the MPI module rather than unchecked external declarations.
      USE mpi
      USE mpiu_I, ONLY: startmpi, mpix_bytes
      IMPLICIT NONE

      character*(MPI_MAX_PROCESSOR_NAME) host
      INTEGER myid, nprocs, lenhost
      INTEGER ierr, ierrtotal
      INTEGER istat(MPI_STATUS_SIZE)



      END MODULE mpi_C
