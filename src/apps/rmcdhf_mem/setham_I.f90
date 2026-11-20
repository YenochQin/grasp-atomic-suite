      module setham_i
      implicit none (type, external)
      interface
      subroutine setham (symmetry_block_index, mpi_rank, mpi_process_count, total_contribution_count, &
         & total_secondary_contribution_count)
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none (type, external)
      integer, intent(in) :: symmetry_block_index
      integer, intent(in) :: mpi_rank
      integer, intent(in) :: mpi_process_count
      integer(int64) :: total_contribution_count(20), total_secondary_contribution_count(20)
      end subroutine setham
      end interface
      end module setham_i
