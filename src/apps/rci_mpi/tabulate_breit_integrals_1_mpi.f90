module tabulate_breit_integrals_1_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: tabulate_breit_integrals_1_mpi
contains

    !***********************************************************************
    !                                                                      *
    subroutine tabulate_breit_integrals_1_mpi (mpi_rank, mpi_process_count, max_orbital_two_j)

        ! Module dependencies
        !   Written by Per Jonsson                Last revision: October 2014  *
        !                                                                      *
        !***********************************************************************
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        use breit_integral_tables, only: breit_integral_indices_1, breit_integral_values_1

        ! Procedure interfaces
        use tabulate_breit_integrals_1_interface, only: tabulate_breit_integrals_1

        ! Module dependencies
        use mpi_utilities, only: allreduce_mpi_real_sum
        use mpi_utilities, only: allreduce_mpi_integer_sum
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: mpi_rank, mpi_process_count
        integer             :: max_orbital_two_j

        ! Local declarations
        integer :: n
        !-----------------------------------------------
        !
        call tabulate_breit_integrals_1 ((mpi_rank), (mpi_process_count), n, max_orbital_two_j)

        ! Gather integrals (and their indeces) from- and send to- all nodes

        call allreduce_mpi_integer_sum (breit_integral_indices_1, n)
        call allreduce_mpi_real_sum (breit_integral_values_1, n)

        return
    end subroutine tabulate_breit_integrals_1_mpi
end module tabulate_breit_integrals_1_mpi_interface
