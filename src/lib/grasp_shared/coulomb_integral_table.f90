module coulomb_integral_table
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: ndcoea, ncoei
    integer, pointer :: indoei(:)
    real(real64), pointer :: valoei(:)
    logical :: frstco
end module coulomb_integral_table
