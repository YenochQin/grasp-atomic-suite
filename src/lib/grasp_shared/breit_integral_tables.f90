module breit_integral_tables
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: ndtpa(6), ntpi(6)
    integer, pointer :: indtp1(:), indtp2(:), indtp3(:), indtp4(:), &
        indtp5(:), indtp6(:)
    real(real64), pointer :: valtp1(:), valtp2(:), valtp3(:), valtp4(:), &
        valtp5(:), valtp6(:)
    logical :: first(6)
end module breit_integral_tables
