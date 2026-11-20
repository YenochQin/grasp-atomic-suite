module radial_quadrature_workspace
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only:  radial_workspace_capacity
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    ! quad reads ta as f(r)*dr/dt through mtp and uses the next three entries as padding.
    ! Other radial routines reuse ta and tb for intermediate arrays between integrations.
    integer :: mtp
    real(real64) :: ta(radial_workspace_capacity), tb(radial_workspace_capacity)
end module radial_quadrature_workspace
