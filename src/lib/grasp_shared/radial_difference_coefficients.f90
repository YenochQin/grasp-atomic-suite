module radial_difference_coefficients
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    ! Partial-panel Newton-Cotes weights are indexed by sample and upper integration node.
    ! The five-point table omits the first sample because the Hartree integrand vanishes there.
    real(real64) :: cnc5c(2:5, 2:4)
    real(real64) :: cnc6c(6, 2:6)
end module radial_difference_coefficients
