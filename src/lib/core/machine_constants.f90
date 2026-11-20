! Set machine precision and exponential limits from LAPACK's floating-point queries.
module machine_constants
    implicit none (type, external)

    ! Module interface
    private
    public :: setmc
contains
    subroutine setmc
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
        use atomic_calculation_state, only: maximum_exp_argument => expmax, minimum_exp_argument => expmin, &
            & machine_precision => precis, maximum_decimal_exponent => tenmax
        use blas_lapack_interfaces, only: dlamch
        implicit none (type, external)

        ! Local declarations
        real(real64) :: machine_limit
        maximum_decimal_exponent = dlamch("L")
        machine_limit = dlamch("O")
        maximum_exp_argument = log(machine_limit)
        machine_limit = dlamch("U")
        minimum_exp_argument = log(machine_limit)
        machine_precision = dlamch("E")
        if (ldbpg(1)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,"From SUBROUTINE SETMC:"/," TENMAX (maximum exponent of 10): ",f5.0,/,&
                &" EXPMAX (maximum exponent of e): ",1p,1d19.12,/," EXPMIN (minimum exponent of e): ",1d19.12,/,&
                &" PRECIS (machine precision): ",1d19.12)') &
                maximum_decimal_exponent, maximum_exp_argument, minimum_exp_argument, machine_precision
        end if
        return
        return
    end subroutine setmc
end module machine_constants
