module tabulate_slater_integrals_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: tabulate_slater_integrals_mpi
contains

    !***********************************************************************
    !                                                                      *
    subroutine tabulate_slater_integrals_mpi(mpi_rank, mpi_process_count, max_orbital_two_j)
        ! Module dependencies
        !
        !   Written by     Xinghong He            Last revision: 12 Jun 1998   *
        !                                                                      *
        !***********************************************************************
        use distributed_slater_integrals, only: local_slater_integral_indices, local_slater_integral_values

        ! Procedure interfaces
        use tabulate_slater_integrals_interface, only: tabulate_slater_integrals

        ! Module dependencies
        use mpi_utilities, only: allreduce_mpi_real_sum
        use mpi_utilities, only: allreduce_mpi_integer_sum
        implicit none (type, external)

        ! Arguments
        integer  :: mpi_rank
        integer  :: mpi_process_count
        integer  :: max_orbital_two_j

        ! Local declarations
        integer      :: n
        !-----------------------------------------------

        call tabulate_slater_integrals (mpi_rank, mpi_process_count, n, max_orbital_two_j)

        ! Gather integrals (and their indeces) from- and send to- all nodes

        call allreduce_mpi_integer_sum (local_slater_integral_indices, n)
        call allreduce_mpi_real_sum (local_slater_integral_values, n)

        return
    end subroutine tabulate_slater_integrals_mpi
end module tabulate_slater_integrals_mpi_interface
