module orbital_lagrange_multipliers
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: nec
    real(real64), pointer :: ecv(:)
    integer, pointer :: iecc(:)
end module orbital_lagrange_multipliers
