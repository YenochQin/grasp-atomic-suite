!***********************************************************************
!                                                                      *
subroutine build_vacuum_polarization
    ! Module dependencies
    !                                                                      *
    !   This  routine controls the setting up of the vacuum polarization   *
    !   potential for the given nuclear charge distribution at each grid   *
    !   point using the analytic functions defined by  L Wayne Fullerton   *
    !   and G A Rinker Jr in Phys Rev A 13 (1976) 1283. The potential is   *
    !   accumulated  in  array TB(I), I = 1, ..., N  which  is in COMMON   *
    !   block /TATB/ .                                                     *
    !                                                                      *
    !   Call(s) to: [RCI92]: VAC2, VAC4.                                   *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
    use nuclear_charge_distribution, only: max_radial_points, zdist
    use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb

    ! Procedure interfaces
    use build_second_order_vacuum_polarization_interface, only: build_second_order_vacuum_polarization
    use build_fourth_order_vacuum_polarization_interface, only: build_fourth_order_vacuum_polarization
    implicit none (type, external)

    ! Local declarations
    integer :: i
    !-----------------------------------------------
    !
    !   Redefine  ZDIST  to be  rho*r*r'
    !
    zdist(:mtp) = zdist(:mtp)*radial_grid(:mtp)*radial_grid_derivative(:mtp)
    !
    !   Second-order vacuum polarisation potential; returned in
    !   array TB
    !
    call build_second_order_vacuum_polarization
    !
    !   Fourth-order vacuum polarization potential; returned in
    !   array TA
    !
    call build_fourth_order_vacuum_polarization
    !
    !   If option 7 is set, use user-defined vacuum polarization
    !   potential
    !
    return
end subroutine build_vacuum_polarization
