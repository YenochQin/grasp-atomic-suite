module initialize_sparse_eigenvectors_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine initialize_sparse_eigenvectors (nmax,csf_count, mpi_rank,mpi_process_count,initial_vector_count,basis, &
            & restart_matrix_unit, &
            & hamiltonian_energy_shift)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: nmax, csf_count, mpi_rank, mpi_process_count, initial_vector_count, restart_matrix_unit
            real(real64), intent(in) :: hamiltonian_energy_shift
            real(real64), intent(inout):: basis(*)
        end subroutine initialize_sparse_eigenvectors
    end interface
end module initialize_sparse_eigenvectors_interface
