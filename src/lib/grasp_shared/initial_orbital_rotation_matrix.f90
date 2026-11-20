module initial_orbital_rotation_matrix
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module constants
    integer, parameter :: nlmax = 20

    ! Module state
    real(real64) :: cici(20*nlmax*nlmax)
end module initial_orbital_rotation_matrix
