module breit_frequency_scale
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: transverse_photon_frequency_scale
end module breit_frequency_scale
