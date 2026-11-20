module mpi_calculation_state
    ! Module dependencies
    use calculation_state_files, only: select_state
    use mpi_state, only: mpi_rank, mpi_bcast, mpi_character, mpi_logical, mpi_comm_world
    implicit none (type, external)

    ! Module interface
    private
    public :: select_mpi_state
contains
    ! Only rank 0 consumes the input record. Every rank receives the resolved name and answer.
    subroutine select_mpi_state(name, question, answer, serial_directory)
        ! Arguments
        character(*), intent(out) :: name
        character(*), intent(in) :: question, serial_directory
        logical, intent(out) :: answer

        ! Local declarations
        integer :: status
        if (mpi_rank == 0) call select_state(name, question, answer, serial_directory)
        call mpi_bcast(name, len(name), mpi_character, 0, mpi_comm_world, status)
        call mpi_bcast(answer, 1, mpi_logical, 0, mpi_comm_world, status)
    end subroutine
end module mpi_calculation_state
