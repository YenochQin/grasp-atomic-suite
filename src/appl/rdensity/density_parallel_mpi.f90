! MPI processes isolate the mutable angular-library module state.
module density_parallel
    use vast_kind_param, only: DOUBLE
    use mpi
    use mpi_C, only: rank => MYID, workers => NPROCS, HOST, LENHOST
    implicit none
    private
    public :: rank, workers, parallel_start, parallel_finish
    public :: broadcast_options, sum_density, density_fail
    integer :: start_count
    interface
        subroutine MPIX_STARTUP(myid, nprocs, host, lenhost, ncount1, &
                progname, description, infiledesc, outfiledesc, quiet_workers)
            integer :: myid, nprocs, lenhost, ncount1
            character(*), intent(in) :: host, progname, description, infiledesc, outfiledesc
            logical, intent(in) :: quiet_workers
        end subroutine
        subroutine MPIX_SHUTDOWN(myid, ncount1, progname)
            integer, intent(in) :: myid, ncount1
            character(*), intent(in) :: progname
        end subroutine
        subroutine GDRSUMMPI_ROOT(x, n)
            import DOUBLE
            integer, intent(in) :: n
            real(DOUBLE), intent(inout) :: x(n)
        end subroutine
    end interface
contains
    subroutine parallel_start
        call MPIX_STARTUP(rank, workers, HOST, LENHOST, start_count, &
            'RDENSITY_MPI', 'MPI density and natural orbitals', &
            'isodata, name.c, name.w, name.(c)m', 'name.nw, name.(c)d', .true.)
    end subroutine

    subroutine parallel_finish
        call MPIX_SHUTDOWN(rank, start_count, 'RDENSITY_MPI')
    end subroutine

    subroutine broadcast_options(name, nci, ndef, ordering, status)
        character(*), intent(inout) :: name
        integer, intent(inout) :: nci, ndef, ordering, status
        integer :: ierr
        call MPI_Bcast(name, len(name), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(nci, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(ndef, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(ordering, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    end subroutine

    subroutine sum_density(matrix)
        real(DOUBLE), contiguous, intent(inout) :: matrix(:,:)
        ! Reduce one state's NW*NW block at a time, using the suite helper.
        call GDRSUMMPI_ROOT(matrix, size(matrix))
    end subroutine

    subroutine density_fail(message)
        character(*), intent(in) :: message
        integer :: ierr
        write(*,*) 'RDENSITY rank ', rank, ': ', message
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        error stop 1
    end subroutine
end module density_parallel
