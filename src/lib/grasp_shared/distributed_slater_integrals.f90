! Shares the Slater integral cache and counts each rank's share of integral evaluations.
module distributed_slater_integrals
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer, pointer :: local_slater_integral_indices(:)
    real(real64), pointer :: local_slater_integral_values(:)

contains

    ! Slater integral tabulation assigns one-based entry i to rank mod(i, mpi_process_count).
    pure integer function rk_local_count(integral_count, mpi_rank, mpi_process_count) result(local_integral_count)
        ! Arguments
        integer, intent(in) :: integral_count, mpi_rank, mpi_process_count
        if (integral_count < 0 .or. mpi_process_count < 1) error stop "RK_LOCAL_COUNT: invalid rank layout"
        if (mpi_rank < 0 .or. mpi_rank >= mpi_process_count) error stop "RK_LOCAL_COUNT: invalid rank layout"
        local_integral_count = integral_count/mpi_process_count
        if (mpi_rank > 0 .and. mpi_rank <= mod(integral_count, mpi_process_count)) &
            local_integral_count = local_integral_count + 1
    end function rk_local_count
end module distributed_slater_integrals
