! Set machine precision and exponential limits from LAPACK's floating-point queries.
module machine_constants
    implicit none (type, external)
    private
    public :: setmc
contains
      subroutine setmc
      use, intrinsic :: iso_fortran_env, only: real64
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use atomic_calculation_state, only: maximum_exp_argument => expmax, minimum_exp_argument => expmin, &
         & machine_precision => precis, maximum_decimal_exponent => tenmax
      use blas_lapack_interfaces, only: dlamch
      implicit none (type, external)
      real(real64) :: machine_limit
      maximum_decimal_exponent = dlamch("L")
      machine_limit = dlamch("O")
      maximum_exp_argument = log(machine_limit)
      machine_limit = dlamch("U")
      minimum_exp_argument = log(machine_limit)
      machine_precision = dlamch("E")
      if (ldbpg(1)) write (99, '(/,"From SUBROUTINE SETMC:"/," TENMAX (maximum exponent of 10): ",f5.0,/,&
         &" EXPMAX (maximum exponent of e): ",1p,1d19.12,/," EXPMIN (minimum exponent of e): ",1d19.12,/,&
         &" PRECIS (machine precision): ",1d19.12)') &
         maximum_decimal_exponent, maximum_exp_argument, minimum_exp_argument, machine_precision
      return
      return
      end subroutine setmc
end module machine_constants
