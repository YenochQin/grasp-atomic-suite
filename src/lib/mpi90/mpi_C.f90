      module mpi_c
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  06:23:52  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
! Use the MPI implementation's explicit Fortran interfaces. Integer handles
! remain the public convention of this legacy module; choice buffers are
! described by the MPI module rather than unchecked external declarations.
      use mpi
      use mpiu_i, only: startmpi, mpix_bytes
      implicit none

      character*(mpi_max_processor_name) host
      integer myid, nprocs, lenhost
      integer ierr, ierrtotal
      integer istat(mpi_status_size)



      end module mpi_c
