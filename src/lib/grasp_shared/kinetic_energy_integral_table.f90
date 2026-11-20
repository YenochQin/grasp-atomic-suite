module kinetic_energy_integral_table
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: ndkea, nkei
    integer, pointer :: indkei(:)
    real(real64), pointer :: valkei(:)
    logical :: frstki
end module kinetic_energy_integral_table
