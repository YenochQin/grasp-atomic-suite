! Combine the radial vacuum-polarization contributions used in RCI.
module build_vacuum_polarization_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: build_vacuum_polarization
contains

    !   This  routine controls the setting up of the vacuum polarization
    !   potential for the given nuclear charge distribution at each grid
    !   point using the analytic functions defined by  L Wayne Fullerton
    !   and G A Rinker Jr in Phys Rev A 13 (1976) 1283. The potential is
    !   accumulated  in  array TB(I), I = 1, ..., N  which  is in COMMON
    !   block /TATB/ .
    !
    ! Written by Farid A Parpia, at Oxford. Last update: 09 Oct 1992
    subroutine build_vacuum_polarization
        ! Module dependencies
        use radial_grid_state, only: radial_grid, radial_grid_derivative
        use nuclear_charge_distribution, only: zdist
        use radial_quadrature_workspace, only: integration_last_index

        ! Procedure interfaces
        use build_second_order_vacuum_polarization_interface, only: build_second_order_vacuum_polarization
        use build_fourth_order_vacuum_polarization_interface, only: build_fourth_order_vacuum_polarization
        implicit none (type, external)

        ! Local declarations
        !
        !   Redefine  ZDIST  to be  rho*r*r'
        zdist(:integration_last_index) = &
            zdist(:integration_last_index)*radial_grid(:integration_last_index)*radial_grid_derivative(:integration_last_index)
        ! Second-order vacuum polarisation potential; returned in
        ! array TB
        call build_second_order_vacuum_polarization
        ! Fourth-order vacuum polarization potential; returned in
        ! array TA
        call build_fourth_order_vacuum_polarization
        ! If option 7 is set, use user-defined vacuum polarization
        ! potential
        return
    end subroutine build_vacuum_polarization
end module build_vacuum_polarization_interface
