module tabulate_breit_integrals_2_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        subroutine tabulate_breit_integrals_2 (mpi_rank, mpi_process_count, nb, max_orbital_two_j)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: mpi_rank, mpi_process_count
            integer             :: nb, max_orbital_two_j
        end subroutine tabulate_breit_integrals_2
    end interface
end module tabulate_breit_integrals_2_interface
