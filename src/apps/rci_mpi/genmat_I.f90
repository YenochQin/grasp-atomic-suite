      module genmat_i
      implicit none (type, external)
      interface
      subroutine genmat (atwinv, symmetry_block_index, mpi_rank, mpi_process_count, elsto, irestart, slf_en)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64) :: atwinv
      integer, intent(in) :: symmetry_block_index
      integer, intent(in) :: mpi_rank
      integer, intent(in) :: mpi_process_count
      real(real64), intent(out) :: elsto
      integer, intent(out) :: irestart
      real(real64) :: slf_en(*)
      end subroutine genmat
      end interface
      end module genmat_i
