module radial_derivative_stencil
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    ! a(i,j) weights sample j for the step-scaled derivative at node i of a 13-point stencil.
    real(real64) :: derivative_weights(13, 13)
end module radial_derivative_stencil
