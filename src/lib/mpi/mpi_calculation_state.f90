module mpi_calculation_state
    ! Module dependencies
    use calculation_state_files, only: read_state_name, read_state_pair, read_state_settings
    use mpi_io_errors, only: abort_input_mpi
    use mpi_state, only: mpi_rank, mpi_bcast, mpi_character, mpi_logical, mpi_comm_world
    implicit none (type, external)

    ! Module interface
    private
    public :: read_mpi_state_name, read_mpi_state_pair, read_mpi_state_settings
contains
    ! Only rank 0 reads the input. A failed validation aborts before broadcasting a name.
    subroutine read_mpi_state_name(name)
        ! Arguments
        character(*), intent(out) :: name

        ! Local declarations
        integer :: status
        if (mpi_rank == 0) then
            call read_state_name(name, status=status)
            if (status /= 0) call abort_input_mpi()
        end if
        call mpi_bcast(name, len(name), mpi_character, 0, mpi_comm_world, status)
    end subroutine

    subroutine read_mpi_state_pair(names)
        ! Arguments
        character(*), intent(out) :: names(2)

        ! Local declarations
        integer :: status
        if (mpi_rank == 0) then
            call read_state_pair(names(1), names(2), status=status)
            if (status /= 0) call abort_input_mpi()
        end if
        call mpi_bcast(names, size(names)*len(names), mpi_character, 0, mpi_comm_world, status)
    end subroutine

    subroutine read_mpi_state_settings(question, answer)
        ! Arguments
        character(*), intent(in) :: question
        logical, intent(out) :: answer

        ! Local declarations
        integer :: status
        if (mpi_rank == 0) call read_state_settings(question, answer)
        call mpi_bcast(answer, 1, mpi_logical, 0, mpi_comm_world, status)
    end subroutine
end module mpi_calculation_state
