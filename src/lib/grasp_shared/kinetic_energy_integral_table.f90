module kinetic_energy_integral_table
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: kinetic_integral_capacity, kinetic_integral_count
    integer, pointer :: kinetic_integral_indices(:)
    real(real64), pointer :: kinetic_integral_values(:)
    logical :: first_kinetic_integral_lookup
end module kinetic_energy_integral_table
