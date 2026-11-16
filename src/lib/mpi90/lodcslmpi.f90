module lodcslmpi_i
    implicit none (type, external)
    private
    public :: lodcslmpi
contains
!***********************************************************************
      subroutine lodcslmpi (nfile, ncore, jblock)

! An MPI container of lodcsh2 which loads CSL list of the current block
! into memory. It forwards the call together with the same set of
! parameters to lodcsh2 and then broadcasts the results to all nodes.
!
! Note: Memories have been allocated/deallocated each block outside.
! This subroutine calls lodcsh2 on node-0 to generate the data for the
! block; and then broadcasts to all other nodes. A new MPI data type
! of 4 byte-long is created to handle 64-bit machines whose MPI
! implementation does not support 4-byte integers. If jblock=-119,
! then ALL blocks will be loaded instead of just one. This is
! implemented in lodcsh2.
!
! Currently used by rcimpivu, mcpmpi, rscfmpivu
!
! Xinghong He 98-08-06
!
!***********************************************************************
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def,   only:  max_orbitals
      use mpi_c
      use orb_c,           only: csf_count, csf_subshell_occupations
      use stat_c,          only: csf_subshell_quantum_codes, csf_coupling_quantum_codes

      use lodcsh2_i, only: lodcsh2
      implicit none (type, external)
      integer       :: nfile, ncore, jblock
      integer :: mpix_int1

!-----------------------------------------------------------------------

      if (mpi_rank == 0) then
         call lodcsh2 ((nfile), (ncore), (jblock))
      end if

! Construct mpi data type for Integer*1 and then broadcast.

!cjb mpix_bytes
!cjb if your compiler or your MPI does not accept type MPIX_INT1
!cjb try standard MPI type MPI_INTEGER1
      call mpix_bytes (1, mpix_int1, mpi_error_code)
!     CALL MPI_Type_create_f90_integer (15, MPIX_INT1, ierr)

      call mpi_bcast                                                   &
               (csf_subshell_occupations(:,:),     max_orbitals*csf_count,mpix_int1,0,mpi_comm_world,mpi_error_code)
!cjb mpix_bytes
!cjb           (IQA(:,:),    NNNW*NCF,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast                                                   &
               (csf_subshell_quantum_codes(:,:,:),3*max_orbitals*csf_count,mpix_int1,0,mpi_comm_world,mpi_error_code)
!cjb mpix_bytes
!cjb           (JQSA(:,:,:),3*NNNW*NCF,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast                                                   &
               (csf_coupling_quantum_codes(:,:),   max_orbitals*csf_count,mpix_int1,0,mpi_comm_world,mpi_error_code)
!cjb mpix_bytes
!cjb           (JCUPA(:,:),   NNNW*NCF,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
      return
      end subroutine lodcslmpi
end module lodcslmpi_i
