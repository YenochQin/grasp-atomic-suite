module coulomb_integral_table
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: coulomb_integral_capacity, coulomb_integral_count
    integer, pointer :: coulomb_integral_indices(:)
    real(real64), pointer :: coulomb_integral_values(:)
    logical :: first_coulomb_integral_lookup
end module coulomb_integral_table
