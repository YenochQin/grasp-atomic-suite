! MPI processes isolate the mutable angular-library module state.
module density_parallel
    use mpi_csf_basis_input, only: load_density_basis => load_csf_basis_mpi
    use, intrinsic :: iso_fortran_env, only: real64
    use mpi, only: mpi_abort, mpi_bcast, mpi_character, mpi_comm_world, mpi_integer
    use mpi_utilities, only: mpix_startup, mpix_shutdown, gdrsummpi_root
    use mpi_state, only: rank => mpi_rank, workers => mpi_process_count, processor_name, processor_name_length
    implicit none (type, external)
    private
    public :: rank, workers, parallel_start, parallel_finish
    public :: broadcast_options, sum_density, density_fail
    public :: load_density_basis
    integer :: start_count
contains
    subroutine parallel_start
        call mpix_startup(rank, workers, processor_name, processor_name_length, start_count, &
            "RDENSITY_MPI", "MPI density and natural orbitals", &
            "isodata, descriptor Parquet/header, name.wfn, name.mix/name.cimix", "name.nw, name.(c)d", .true.)
    end subroutine parallel_start

    subroutine parallel_finish
        call mpix_shutdown(rank, start_count, "RDENSITY_MPI")
    end subroutine parallel_finish

    subroutine broadcast_options(name, nci, ndef, ordering, status)
        character(*), intent(inout) :: name
        integer, intent(inout) :: nci, ndef, ordering, status
        integer :: ierr
        call mpi_bcast(name, len(name), mpi_character, 0, mpi_comm_world, ierr)
        call mpi_bcast(nci, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(ndef, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(ordering, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(status, 1, mpi_integer, 0, mpi_comm_world, ierr)
    end subroutine broadcast_options

    subroutine sum_density(matrix)
        real(real64), contiguous, intent(inout) :: matrix(:,:)
        ! Reduce one state's NW*NW block at a time, using the suite helper.
        call gdrsummpi_root(matrix, size(matrix))
    end subroutine sum_density

    subroutine density_fail(message)
        character(*), intent(in) :: message
        integer :: ierr
        write(*,*) "RDENSITY rank ", rank, ": ", message
        call mpi_abort(mpi_comm_world, 1, ierr)
        error stop 1
    end subroutine density_fail
end module density_parallel
