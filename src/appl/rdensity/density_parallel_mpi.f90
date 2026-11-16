! MPI processes isolate the mutable angular-library module state.
module density_parallel
    use, intrinsic :: iso_fortran_env, only: real64
    use mpi
    use mpiu_i, only: mpix_startup, mpix_shutdown, gdrsummpi_root
    use mpi_c, only: rank => mpi_rank, workers => mpi_process_count, processor_name, processor_name_length
    implicit none (type, external)
    private
    public :: rank, workers, parallel_start, parallel_finish
    public :: broadcast_options, sum_density, density_fail
    integer :: start_count
contains
    subroutine parallel_start
        call mpix_startup(rank, workers, processor_name, processor_name_length, start_count, &
            'RDENSITY_MPI', 'MPI density and natural orbitals', &
            'isodata, name.c, name.w, name.(c)m', 'name.nw, name.(c)d', .true.)
    end subroutine

    subroutine parallel_finish
        call mpix_shutdown(rank, start_count, 'RDENSITY_MPI')
    end subroutine

    subroutine broadcast_options(name, nci, ndef, ordering, status)
        character(*), intent(inout) :: name
        integer, intent(inout) :: nci, ndef, ordering, status
        integer :: ierr
        call mpi_bcast(name, len(name), mpi_character, 0, mpi_comm_world, ierr)
        call mpi_bcast(nci, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(ndef, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(ordering, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(status, 1, mpi_integer, 0, mpi_comm_world, ierr)
    end subroutine

    subroutine sum_density(matrix)
        real(real64), contiguous, intent(inout) :: matrix(:,:)
        ! Reduce one state's NW*NW block at a time, using the suite helper.
        call gdrsummpi_root(matrix, size(matrix))
    end subroutine

    subroutine density_fail(message)
        character(*), intent(in) :: message
        integer :: ierr
        write(*,*) 'RDENSITY rank ', rank, ': ', message
        call mpi_abort(mpi_comm_world, 1, ierr)
        error stop 1
    end subroutine
end module density_parallel
