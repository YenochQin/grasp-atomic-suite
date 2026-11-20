!***********************************************************************
!                                                                      *
subroutine tabulate_slater_integrals_mpi(mpi_rank, mpi_process_count, max_orbital_two_j)
    ! Module dependencies
    !
    !   Written by     Xinghong He            Last revision: 12 Jun 1998   *
    !                                                                      *
    !***********************************************************************
    use distributed_slater_integrals, only: indteirk, valteirk, rk_local_count

    ! Procedure interfaces
    use tabulate_slater_integrals_interface, only: tabulate_slater_integrals

    ! Module dependencies
    use mpi_utilities, only: gdsummpi
    use mpi_utilities, only: gisummpi
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

    call gisummpi (indteirk, n)
    call gdsummpi (valteirk, n)

    return
end subroutine tabulate_slater_integrals_mpi
