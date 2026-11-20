! Provides the suite's standard input, output, and error unit numbers.
module standard_io_units
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: input_unit, output_unit, error_unit
    implicit none (type, external)

    ! Module state
    integer :: stdin_unit=input_unit, stdout_unit=output_unit, stderr_unit=error_unit
end module standard_io_units
