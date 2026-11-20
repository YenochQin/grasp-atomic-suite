! Populate the table of log((index-1)!) used by angular-momentum routines.
module log_factorial_table
    implicit none (type, external)

    ! Module interface
    private
    public :: initialize_log_factorial_table
contains
    subroutine initialize_log_factorial_table
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use log_factorial_cache, only: log_factorials, log_factorial_table_size
        implicit none (type, external)

        ! Local declarations
        integer :: factorial_table_index
        real(real64) :: factorial_multiplier
        log_factorials(1) = 1.0_real64
        log_factorials(2) = 1.0_real64
        factorial_multiplier = 2.0_real64
        do factorial_table_index = 3, 30
            log_factorials(factorial_table_index) = log_factorials(factorial_table_index-1)*factorial_multiplier
            factorial_multiplier = factorial_multiplier + 1.0_real64
        end do
        do factorial_table_index = 1, 30
            log_factorials(factorial_table_index) = log(log_factorials(factorial_table_index))
        end do
        factorial_multiplier = 3.0e1_real64
        do factorial_table_index = 31, log_factorial_table_size
            log_factorials(factorial_table_index) = log_factorials(factorial_table_index-1) + log(factorial_multiplier)
            factorial_multiplier = factorial_multiplier + 1.0_real64
        end do
        return
    end subroutine initialize_log_factorial_table
end module log_factorial_table
