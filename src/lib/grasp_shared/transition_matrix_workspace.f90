module transition_matrix_workspace
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: nsdim
    real(real64), pointer :: hb1(:), hb2(:), hc1(:), hc2(:), hm1(:), hm2(:)
    integer, pointer :: jja(:), jjb(:)
    integer :: lk, kk
    integer :: ntdim
    real(real64), pointer :: xsldr(:), totc(:), totb(:)
    integer, pointer :: isldr(:), isldr1(:)
    integer :: nint, ninteg
    integer, pointer :: nptr(:), lab(:)
    integer :: nkp
    integer, pointer :: kp(:)
    logical :: ltc(10)
    integer :: initial_mixing_unit, final_mixing_unit
end module transition_matrix_workspace
