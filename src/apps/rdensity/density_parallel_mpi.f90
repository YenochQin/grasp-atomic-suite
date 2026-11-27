! Provide MPI option broadcasts, density reductions, and error handling for rdensity_mpi.
! MPI processes isolate the mutable angular-library module state.
module density_parallel
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use mpi_f08, only: mpi_bcast, mpi_character, mpi_comm_world, mpi_count_kind, mpi_integer, &
        mpi_double_precision, mpi_in_place, mpi_reduce, mpi_success, mpi_sum
    use mpi_element_counts, only: checked_mpi_count
    use mpi_io_errors, only: density_abort => abort_input_mpi
    use mpi_utilities, only: start_mpi_application, finish_mpi_application
    use mpi_state, only: mpi_rank, mpi_process_count, processor_name, processor_name_length
    implicit none (type, external)

    ! Module interface
    private
    public :: mpi_rank, mpi_process_count, parallel_start, parallel_finish
    public :: broadcast_options, sum_density, density_fail, density_abort
    public :: load_density_basis

    ! Module state
    integer(int64) :: start_count
contains
    subroutine load_density_basis(descriptor_path, core_subshell_count)
        ! Module dependencies
        use mpi_csf_basis_input, only: load_csf_basis_mpi

        ! Arguments
        character(len=*), intent(in) :: descriptor_path
        integer, intent(out) :: core_subshell_count

        call load_csf_basis_mpi(descriptor_path, core_subshell_count)
    end subroutine load_density_basis

    subroutine parallel_start
        call start_mpi_application(mpi_rank, mpi_process_count, processor_name, processor_name_length, start_count, &
            "RDENSITY_MPI", "MPI density and natural orbitals", &
            "isodata, descriptor Parquet/header, name.wfn, name.mix/name.cimix", &
            "name.natural.wfn or name.ci.natural.wfn, name.(c)d", .true.)
    end subroutine parallel_start

    subroutine parallel_finish
        call finish_mpi_application(mpi_rank, start_count, "RDENSITY_MPI")
    end subroutine parallel_finish

    subroutine broadcast_options(name, mixing_file_kind, nondefault_options_flag, ordering, status)
        ! Arguments
        character(*), intent(inout) :: name
        integer, intent(inout) :: mixing_file_kind, nondefault_options_flag, ordering, status

        ! Local declarations
        integer :: mpi_status
        call mpi_bcast(name, len(name), mpi_character, 0, mpi_comm_world, mpi_status)
        call mpi_bcast(mixing_file_kind, 1, mpi_integer, 0, mpi_comm_world, mpi_status)
        call mpi_bcast(nondefault_options_flag, 1, mpi_integer, 0, mpi_comm_world, mpi_status)
        call mpi_bcast(ordering, 1, mpi_integer, 0, mpi_comm_world, mpi_status)
        call mpi_bcast(status, 1, mpi_integer, 0, mpi_comm_world, mpi_status)
    end subroutine broadcast_options

    subroutine sum_density(matrix)
        ! Arguments
        real(real64), contiguous, intent(inout) :: matrix(:,:)

        ! Local declarations
        integer :: element_count, mpi_status
        real(real64) :: unused_receive_buffer(1)

        element_count = checked_mpi_count(size(matrix,kind=mpi_count_kind), "RDENSITY matrix reduction")
        if (mpi_rank == 0) then
            call mpi_reduce(mpi_in_place, matrix, element_count, mpi_double_precision, mpi_sum, 0, mpi_comm_world, mpi_status)
        else
            ! MPI ignores the receive buffer on workers, including ranks with no CSF columns.
            call mpi_reduce(matrix, unused_receive_buffer, element_count, mpi_double_precision, mpi_sum, 0, &
                mpi_comm_world, mpi_status)
        end if
        if (mpi_status /= mpi_success) call density_fail("Density matrix reduction failed")
    end subroutine sum_density

    subroutine density_fail(message)
        ! Arguments
        character(*), intent(in) :: message

        write(*,*) "RDENSITY rank ", mpi_rank, ": ", message
        call density_abort()
        error stop 1
    end subroutine density_fail
end module density_parallel
