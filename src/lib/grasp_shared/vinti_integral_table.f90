module vinti_integral_table
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: ndvin, nvinti
    integer, pointer :: indtei(:)
    real(real64), pointer :: valtei(:)
    logical :: frstvi
end module vinti_integral_table
