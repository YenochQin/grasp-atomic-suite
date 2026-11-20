module rci_prerun_settings
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: iprerun, ncsfpre
    real(real64) :: coeffcut1, coeffcut2
end module rci_prerun_settings
