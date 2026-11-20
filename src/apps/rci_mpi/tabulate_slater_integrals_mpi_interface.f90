module tabulate_slater_integrals_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine tabulate_slater_integrals_mpi (mpi_rank, mpi_process_count, max_orbital_two_j)
            implicit none (type, external)

            ! Arguments
            integer :: mpi_rank
            integer :: mpi_process_count
            integer :: max_orbital_two_j
        end subroutine tabulate_slater_integrals_mpi
    end interface
end module tabulate_slater_integrals_mpi_interface
