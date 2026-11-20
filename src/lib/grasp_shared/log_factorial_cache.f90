! Supplies the logarithmic factorial table used by the angular-momentum algebra.
module log_factorial_cache
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module constants
    integer, parameter :: log_factorial_table_size = 500

    ! Module state
    ! factt initializes gam(i) = log((i-1)!), so the first slot represents 0!.
    real(real64) :: log_factorials(log_factorial_table_size)
end module log_factorial_cache
