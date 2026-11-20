      module iniestsd_i
      implicit none (type, external)
      interface
      subroutine iniestsd (nmax,csf_count, mpi_rank,mpi_process_count,initial_vector_count,basis,restart_matrix_unit, &
         & hamiltonian_energy_shift)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: nmax, csf_count, mpi_rank, mpi_process_count, initial_vector_count, restart_matrix_unit
      real(real64), intent(in) :: hamiltonian_energy_shift
      real(real64), intent(inout):: basis(*)
      end subroutine iniestsd
      end interface
      end module iniestsd_i
