module radial_wavefunction_derivatives
    implicit none (type, external)
    private
    public :: differentiate_radial_wavefunctions
contains
subroutine differentiate_radial_wavefunctions(j)
    use radial_grid_state, only: radial_point_count
    use radial_quadrature_workspace, only: radial_workspace_a, radial_workspace_b
    implicit none (type, external)
    integer :: j
    radial_workspace_a(:radial_point_count) = 0
    radial_workspace_b(:radial_point_count) = 0
end subroutine differentiate_radial_wavefunctions
end module radial_wavefunction_derivatives

module radial_quadrature
    implicit none (type, external)
    private
    public :: integrate_radial_grid
contains

subroutine integrate_radial_grid(result)
    use, intrinsic :: iso_fortran_env, only: real64
    use radial_quadrature_workspace, only: radial_workspace_a, integration_last_index
    implicit none (type, external)
    real(real64), intent(out) :: result
    integer :: i
    result = 0
    do i = 2, integration_last_index
        result = result+real(i,real64)*radial_workspace_a(i)
    end do
end subroutine integrate_radial_grid
end module radial_quadrature
