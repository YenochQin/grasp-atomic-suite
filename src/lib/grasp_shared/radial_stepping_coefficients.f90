module radial_stepping_coefficients
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: starting_step_weights(6)
end module radial_stepping_coefficients
