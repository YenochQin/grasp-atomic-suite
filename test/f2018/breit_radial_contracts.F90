program breit_radial_contracts
    use, intrinsic :: iso_fortran_env, only: real64
    use breit_radial_integral_coefficients, only: calculate_breit_radial_integral_coefficients
#ifndef CHECK_MPI
    use log_factorial_table, only: initialize_log_factorial_table
#endif
#ifdef CHECK_MPI
    use breit_radial_integral_coefficients, only: set_breit_abort_handler
    use mpi_io_errors, only: abort_input_mpi
    use mpi_state, only: mpi_init, mpi_comm_rank, mpi_comm_world, mpi_barrier, mpi_finalize
#endif
    implicit none (type, external)
    character(16) :: scenario
    integer, parameter :: orbitals(4) = [1,2,1,2], twice_kappa(4) = [-2,-2,-2,-2]
#ifndef CHECK_MPI
    real(real64), parameter :: tolerance = 1.0e-13_real64
#endif
    real(real64) :: factors(12)
#ifdef CHECK_MPI
    integer :: rank, status
#endif

    call get_command_argument(1, scenario)
#ifdef CHECK_MPI
    call mpi_init(status)
    call mpi_comm_rank(mpi_comm_world,rank,status)
    call set_breit_abort_handler(abort_input_mpi)
    if (rank == 1) call invalid_input()
    call mpi_barrier(mpi_comm_world,status)
    call mpi_finalize(status)
#else
    if (scenario /= "valid") call invalid_input()
    call initialize_log_factorial_table()
    call calculate_breit_radial_integral_coefficients(factors,orbitals,twice_kappa,1,1,3,1)
    if (abs(factors(1)-16.0_real64/3) > tolerance .or. any(abs(factors(2:)) > tolerance)) &
        error stop "Direct Breit magnetic coefficient"
    call calculate_breit_radial_integral_coefficients(factors,orbitals,twice_kappa,1,1,4,2)
    if (any(abs(factors(:3)-[4.0_real64,8.0_real64,4.0_real64]/3) > tolerance) .or. &
        any(abs(factors(4:)) > tolerance)) error stop "Exchange Breit magnetic coefficients"
    call calculate_breit_radial_integral_coefficients(factors,orbitals,[-2,-4,-2,-4],1,1,4,2)
    if (any(abs(factors(:3)+[6.0_real64,12.0_real64,6.0_real64]) > tolerance) .or. &
        any(abs(factors(4:)) > tolerance)) error stop "Exchange Breit magnetic phase"
    call calculate_breit_radial_integral_coefficients(factors,orbitals,twice_kappa,2,2,3,1)
    if (any(abs(factors) > tolerance)) error stop "Even tensor rank selection"
#endif
contains
    subroutine invalid_input()
        if (scenario == "rank-zero") then
            call calculate_breit_radial_integral_coefficients(factors,orbitals,twice_kappa,0,0,1,1)
        else
            call calculate_breit_radial_integral_coefficients(factors,[1,2,2,1],twice_kappa,1,1,3,1)
        end if
        error stop "Accepted invalid Breit input"
    end subroutine invalid_input
end program breit_radial_contracts
