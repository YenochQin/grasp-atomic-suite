! Stores the number and values of parameters for the nuclear charge model.
module nuclear_model_parameters_state
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: nuclear_model_parameter_count
    real(real64) :: nuclear_model_parameters(2)
end module nuclear_model_parameters_state
