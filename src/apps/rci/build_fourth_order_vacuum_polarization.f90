!***********************************************************************
!                                                                      *
subroutine build_fourth_order_vacuum_polarization
    ! Module dependencies
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    !                                                                      *
    !   This routine sets up the fourth-order vacuum polarization poten-   *
    !   tial using equations (11) and (12) of L Wayne Fullerton and  G A   *
    !   Rinker, Jr,  Phys  Rev  A 13 (1976) 1283-1287. The potential  is   *
    !   accumulated in array  TC(I), I = 1, ..., N. It is transferred to   *
    !   array TA in COMMON block TATB.                                     *
    !                                                                      *
    !   Call(s) to: [LIB92]: QUAD.                                         *
    !               [RCI92]: FUNL.                                         *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 15 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_radial_points
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use atomic_calculation_state, only: max_radial_points, nuclear_charge, pi, precis, vacuum_speed_of_light_au
    use radial_grid_state, only: radial_grid, radial_point_count
    use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
    use nuclear_charge_distribution, only: zdist
    use radial_quadrature_workspace, only: mtp, ta, tb

    ! Procedure interfaces
    use fullerton_rinker_l_function_interface, only: fullerton_rinker_l_function

    ! Module dependencies
    use radial_quadrature, only: quad
    implicit none (type, external)

    ! Local declarations
    integer :: i, k, nb2, row_count, ii, ii1, ii2
    real(real64) :: tc(max_radial_points)
    real(real64) :: epsi, twocv, factor, ri, x, tci, rk, xk, xi, xm, xp
    !-----------------------------------------------
    !
    !   Overall initialization
    !
    epsi = precis*precis
    twocv = vacuum_speed_of_light_au + vacuum_speed_of_light_au
    !
    !   Potential for point nucleus: equation (12)
    !
    factor = -nuclear_charge/(pi*vacuum_speed_of_light_au)**2
    !
    tc(1) = 0.0_real64
    !
    i = 1
    do
        i = i + 1
        ri = radial_grid(i)
        x = twocv*ri
        tci = (factor/ri)*fullerton_rinker_l_function(x,1)
        if (abs(tci) >= epsi) then
            tc(i) = tci
            if (i >= radial_point_count) exit
        else
            tc(i:radial_point_count) = 0.0_real64
            exit
        end if
    end do
    !
    !   Potential for finite nucleus: equation (11)
    !
    if (nuclear_model_parameter_count == 2) then
        !
        factor = -1.0_real64/(pi*vacuum_speed_of_light_au**3)
        !
        tc(1) = 0.0_real64
        !
        k = 1
        do
            k = k + 1
            !
            rk = radial_grid(k)
            xk = twocv*rk
            ta(1) = 0.0_real64
            !
            do i = 2, mtp
                xi = twocv*radial_grid(i)
                xm = abs(xk - xi)
                xp = xk + xi
                ta(i) = (fullerton_rinker_l_function(xm,0) - fullerton_rinker_l_function(xp,0))*zdist(i)
            end do
            !
            call quad (x)
            !
            x = x*factor/rk
            !
            !   Get out of the loop if the asymptotic region has been reached
            !
            if (abs(x) >= epsi) then
                if (abs((tc(k)-x)/x) > 1.0e-3_real64) then
                    tc(k) = x
                    if (k < radial_point_count) cycle
                end if
            end if
            exit
        end do
        !
    end if
    !
    if (ldbpr(8)) then
        call ensure_debug_output()
        write (angular_debug_unit, '(/,/,/," ++++++++++ VAC4 ++++++++++"/,/,2(" -------- r -------- ----- VV2 (r) -----",&
            &" ----- VV4 (r) -----"))')
        nb2 = radial_point_count/2
        if (2*nb2 == radial_point_count) then
            row_count = nb2
        else
            row_count = nb2 + 1
        end if
        do ii = 1, row_count
            ii1 = ii
            ii2 = ii1 + row_count
            if (ii2 <= radial_point_count) then
                call ensure_debug_output()
                write (angular_debug_unit, '(1p,6(1x,1d19.12))') radial_grid(ii1), tb(ii1), tc(ii1), radial_grid(ii2), tb(ii2), tc(&
                    ii2)
            else if (ii1 <= radial_point_count) then
                call ensure_debug_output()
                write (angular_debug_unit, '(1p,6(1x,1d19.12))') radial_grid(ii1), tb(ii1), tc(ii1)
            end if
        end do
    end if
    !
    !   Generate total vacuum-polarization potential
    !
    tb(:radial_point_count) = tc(:radial_point_count) + tb(:radial_point_count)
    !
    return
    !
    return
    !
end subroutine build_fourth_order_vacuum_polarization
