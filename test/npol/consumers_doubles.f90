module radial_wavefunction_derivatives
    implicit none (type, external)
    private
    public :: dpbdt
contains
subroutine DPBDT(j)
    use radial_grid_state, only: N => radial_point_count
    use radial_quadrature_workspace, only: TA, TB
    implicit none (type, external)
    integer :: j
    TA(:N) = 0
    TB(:N) = 0
end subroutine DPBDT
end module radial_wavefunction_derivatives

module radial_quadrature
    implicit none (type, external)
    private
    public :: quad
contains

subroutine QUAD(result)
    use, intrinsic :: iso_fortran_env, only: real64
    use radial_quadrature_workspace, only: TA, MTP
    implicit none (type, external)
    real(real64), intent(out) :: result
    integer :: i
    result = 0
    do i = 2, MTP
        result = result+real(i,real64)*TA(i)
    end do
end subroutine QUAD
end module radial_quadrature
