module orbital_lagrange_multipliers
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: lagrange_multiplier_count
    real(real64), pointer :: lagrange_multiplier_values(:)
    integer, pointer :: lagrange_orbital_pair_labels(:)
end module orbital_lagrange_multipliers
