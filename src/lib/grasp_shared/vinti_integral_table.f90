module vinti_integral_table
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: vinti_integral_capacity, vinti_integral_count
    integer, pointer :: vinti_integral_indices(:)
    real(real64), pointer :: vinti_integral_values(:)
    logical :: first_vinti_lookup
end module vinti_integral_table
