module build_second_order_vacuum_polarization_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: build_second_order_vacuum_polarization
contains

    !***********************************************************************
    !                                                                      *
    subroutine build_second_order_vacuum_polarization
        ! Module dependencies
        !                                                                      *
        !   This routine sets up the second-order vacuum polarization poten-   *
        !   tial using  equations (1) and (4) of  L Wayne Fullerton and  G A   *
        !   Rinker, Jr,  Phys Rev A  13 (1976) 1283-1287.  The  potential is   *
        !   accumulated  in  array  TB(I), I = 1, ..., N  which is in COMMON   *
        !   block /TATB/ .                                                     *
        !                                                                      *
        !   Call(s) to: [LIB92]: QUAD.                                         *
        !               [RCI92]: FUNK.                                         *
        !                                                                      *
        !   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: nuclear_charge, pi, machine_precision, vacuum_speed_of_light_au
        use radial_grid_state, only: radial_point_count, radial_grid
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count
        use nuclear_charge_distribution, only: zdist
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a, radial_workspace_b

        ! Procedure interfaces
        use fullerton_rinker_k_function_interface, only: fullerton_rinker_k_function

        ! Module dependencies
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Local declarations
        integer :: i, k
        real(real64) :: epsi, twocv, factor, ri, x, tbi, rk, xk, xi, xm, xp
        !-----------------------------------------------
        !
        !   Overall initialization
        !
        epsi = machine_precision*machine_precision
        twocv = vacuum_speed_of_light_au + vacuum_speed_of_light_au
        !
        !   Potential for a point nucleus: equation (1)
        !   (this is also the asymptotoc form for a finite nucleus)
        !
        factor = -(2.0_real64*nuclear_charge)/(3.0_real64*pi*vacuum_speed_of_light_au)
        !
        radial_workspace_b(1) = 0.0_real64
        !
        i = 1
        do
            i = i + 1
            !
            ri = radial_grid(i)
            x = twocv*ri
            tbi = (factor/ri)*fullerton_rinker_k_function(x,1)
            !
            if (abs(tbi) >= epsi) then
                radial_workspace_b(i) = tbi
                if (i >= radial_point_count) exit
            else
                radial_workspace_b(i:radial_point_count) = 0.0_real64
                exit
            end if
        end do
        !
        !   Potential for a finite nucleus: equation (4)
        !
        if (nuclear_model_parameter_count == 2) then
            !
            factor = -2.0_real64/(3.0_real64*vacuum_speed_of_light_au**2)
            !
            !   Set up integrand
            !
            radial_workspace_b(1) = 0.0_real64
            !
            k = 1
            do
                k = k + 1
                !
                rk = radial_grid(k)
                xk = twocv*rk
                !
                radial_workspace_a(1) = 0.0_real64
                do i = 2, integration_last_index
                    xi = twocv*radial_grid(i)
                    xm = abs(xk - xi)
                    xp = xk + xi
                    radial_workspace_a(i) = (fullerton_rinker_k_function(xm,0) - fullerton_rinker_k_function(xp,0))*zdist(i)
                end do
                !
                call integrate_radial_grid (x)
                !
                x = x*factor/rk
                !
                !   Get out of loop if the asymptotic value has been attained
                !
                if (abs(x) >= epsi) then
                    if (abs((x - radial_workspace_b(k))/x) > 1.0e-5_real64) then
                        radial_workspace_b(k) = x
                        if (k < radial_point_count) cycle
                    end if
                end if
                exit
            end do
            !
        end if
        !
        return
    end subroutine build_second_order_vacuum_polarization
end module build_second_order_vacuum_polarization_interface
