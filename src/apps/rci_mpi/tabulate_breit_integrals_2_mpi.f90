!***********************************************************************
!                                                                      *
subroutine tabulate_breit_integrals_2_mpi (mpi_rank, mpi_process_count, max_orbital_two_j)

    ! Module dependencies
    !   Written by Per Jonsson                Last revision: October 2014  *
    !                                                                      *
    !***********************************************************************
    !...Translated by Charlotte Froese Fischer
    !                       Gediminas Gaigalas  10/05/17
    use breit_integral_tables, only: ndtpa, ntpi, indtp1, indtp2, indtp3, indtp4, indtp5, indtp6, &
        valtp1, valtp2, valtp3, valtp4, valtp5, valtp6, first

    ! Procedure interfaces
    use tabulate_breit_integrals_2_interface, only: tabulate_breit_integrals_2

    ! Module dependencies
    use mpi_utilities, only: gdsummpi
    use mpi_utilities, only: gisummpi
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: mpi_rank, mpi_process_count
    integer             :: max_orbital_two_j

    ! Local declarations
    integer :: n
    !-----------------------------------------------
    !
    call tabulate_breit_integrals_2 ((mpi_rank), (mpi_process_count), n, max_orbital_two_j)

    ! Gather integrals (and their indeces) from- and send to- all nodes

    call gisummpi (indtp2, n)
    call gdsummpi (valtp2, n)

    return
end subroutine tabulate_breit_integrals_2_mpi
