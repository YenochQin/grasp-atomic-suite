module radial_quadrature_workspace
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only:  radial_workspace_capacity
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    ! Radial quadrature reads f(r)*dr/dt from radial_workspace_a through integration_last_index.
    ! It uses the next three entries as padding. Other radial routines reuse both arrays for intermediate terms.
    integer :: integration_last_index
    real(real64) :: radial_workspace_a(radial_workspace_capacity), radial_workspace_b(radial_workspace_capacity)
end module radial_quadrature_workspace
