module tabulate_slater_integrals_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine tabulate_slater_integrals (mpi_rank, mpi_process_count, n, max_orbital_two_j)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: mpi_rank, mpi_process_count
            integer, intent(out) :: n, max_orbital_two_j
        end subroutine tabulate_slater_integrals
    end interface
end module tabulate_slater_integrals_interface
