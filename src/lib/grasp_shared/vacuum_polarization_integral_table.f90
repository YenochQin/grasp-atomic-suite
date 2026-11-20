module vacuum_polarization_integral_table
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: ndvpa, nvpi
    integer, pointer :: indvpi(:)
    real(real64), pointer :: valvpi(:)
    logical :: frstvp
end module vacuum_polarization_integral_table
