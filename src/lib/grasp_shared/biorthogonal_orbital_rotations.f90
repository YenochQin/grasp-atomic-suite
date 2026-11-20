module biorthogonal_orbital_rotations
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module constants
    integer, parameter :: nlmax = 20

    ! Module state
    real(real64) :: cfci(20*nlmax*nlmax)
    real(real64) :: cici(20*nlmax*nlmax)
end module biorthogonal_orbital_rotations
