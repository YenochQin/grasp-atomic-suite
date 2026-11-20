module angular_coefficient_buffer_state
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: angular_buffer_capacity, angular_coefficient_count
    ! Labels occupy the first dimension. Column i belongs to angular_coefficients(i).
    integer, pointer :: angular_coefficient_labels(:, :)
    real(real64), pointer :: angular_coefficients(:)
end module angular_coefficient_buffer_state
