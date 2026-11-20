      module setham_i
      implicit none (type, external)
      interface
      subroutine setham (mpi_rank, mpi_process_count, symmetry_block_index, elsto,icstrt, nelmntt,  &
                         atwinv,slf_en)
      use, intrinsic :: iso_fortran_env, only: real64, int64
      implicit none (type, external)
      integer(int64)       :: nelmntt
      integer             :: symmetry_block_index, icstrt
      integer, intent(in) :: mpi_rank, mpi_process_count
      real(real64) :: elsto, atwinv
      real(real64) :: slf_en(*)
      end subroutine setham
      end interface
      end module setham_i
