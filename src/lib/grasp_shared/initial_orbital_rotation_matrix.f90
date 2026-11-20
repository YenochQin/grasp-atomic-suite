module initial_orbital_rotation_matrix
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module constants
    integer, parameter :: rotation_group_capacity = 20

    ! Module state
    real(real64) :: initial_orbital_rotation_values(20*rotation_group_capacity*rotation_group_capacity)
end module initial_orbital_rotation_matrix
