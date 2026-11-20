      module hmoutmpi_i
      implicit none (type, external)
      interface
      subroutine hmoutmpi (mpi_rank, mpi_process_count, block_csf_count)
      implicit none (type, external)
      integer :: mpi_rank
!VAST...Dummy argument MYID is not referenced in this routine.
      integer :: mpi_process_count
!VAST...Dummy argument NPROCS is not referenced in this routine.
      integer, intent(in) :: block_csf_count
      end subroutine hmoutmpi
      end interface
      end module hmoutmpi_i
