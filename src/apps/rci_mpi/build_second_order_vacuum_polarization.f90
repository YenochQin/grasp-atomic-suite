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
    use atomic_calculation_state, only: max_radial_points, nuclear_charge, pi, precis, vacuum_speed_of_light_au
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
    use nuclear_charge_distribution, only: max_radial_points, zdist
    use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb

    ! Procedure interfaces
    use fullerton_rinker_k_function_interface, only: fullerton_rinker_k_function

    ! Module dependencies
    use radial_quadrature, only: quad
    implicit none (type, external)

    ! Local declarations
    integer :: i, k
    real(real64) :: epsi, twocv, factor, ri, x, tbi, rk, xk, xi, xm, xp
    !-----------------------------------------------
    !
    !   Overall initialization
    !
    epsi = precis*precis
    twocv = vacuum_speed_of_light_au + vacuum_speed_of_light_au
    !
    !   Potential for a point nucleus: equation (1)
    !   (this is also the asymptotoc form for a finite nucleus)
    !
    factor = -(2.0_real64*nuclear_charge)/(3.0_real64*pi*vacuum_speed_of_light_au)
    !
    tb(1) = 0.0_real64
    !
    i = 1
    do
        i = i + 1
        !
        ri = radial_grid(i)
        x = twocv*ri
        tbi = (factor/ri)*fullerton_rinker_k_function(x,1)
        !
        if (dabs(tbi) >= epsi) then
            tb(i) = tbi
            if (i >= radial_point_count) exit
        else
            tb(i:radial_point_count) = 0.0_real64
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
        tb(1) = 0.0_real64
        !
        k = 1
        do
            k = k + 1
            !
            rk = radial_grid(k)
            xk = twocv*rk
            !
            ta(1) = 0.0_real64
            do i = 2, mtp
                xi = twocv*radial_grid(i)
                xm = dabs(xk - xi)
                xp = xk + xi
                ta(i) = (fullerton_rinker_k_function(xm,0) - fullerton_rinker_k_function(xp,0))*zdist(i)
            end do
            !
            call quad (x)
            !
            x = x*factor/rk
            !
            !   Get out of loop if the asymptotic value has been attained
            !
            if (dabs(x) >= epsi) then
                if (dabs((x - tb(k))/x) > 1.0e-5_real64) then
                    tb(k) = x
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
