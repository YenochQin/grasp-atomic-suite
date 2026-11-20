module ci_matrix_assembly_state
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: ntpitmp(6)
    integer :: ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
        nvpitmp, nkeitmp, nvintitmp, nelmnttmp, ncftmp
    real(real64) :: cutofftmp
end module ci_matrix_assembly_state
