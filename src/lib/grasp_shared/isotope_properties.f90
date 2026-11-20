module isotope_properties
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: nuclear_spin,nuclear_magnetic_dipole_moment,nuclear_electric_quadrupole_moment
    real(real64) :: hfsi, hfsd, hfsq
    real(real64) :: smsi, smsd, smsq
end module isotope_properties
