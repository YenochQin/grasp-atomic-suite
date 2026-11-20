module isotope_properties
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: sqn,dmomnm,qmomb
    real(real64) :: hfsi, hfsd, hfsq
    real(real64) :: smsi, smsd, smsq
end module isotope_properties
