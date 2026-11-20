module vacuum_polarization_integral_table
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: vacuum_polarization_integral_capacity, vacuum_polarization_integral_count
    integer, pointer :: vacuum_polarization_integral_indices(:)
    real(real64), pointer :: vacuum_polarization_integral_values(:)
    logical :: first_vacuum_polarization_lookup
end module vacuum_polarization_integral_table
