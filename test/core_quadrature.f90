program core_quadrature
    implicit none (type, external)

    call setup
    call runtest

contains

    subroutine runtest
        use grasp_parameters, only: nnn1, radial_workspace_capacity
        use radial_grid_state, only: h => transformed_grid_step, hp => asymptotic_grid_step, n => radial_point_count, &
           & nnn1 => radial_workspace_capacity, r => radial_grid, rnt => radial_grid_scale, rp => radial_grid_derivative, &
           & rpor => radial_derivative_over_radius
        use radial_quadrature_workspace, only: nnn1 => radial_workspace_capacity, mtp, ta
        use radial_quadrature, only: quad
        use, intrinsic :: iso_fortran_env, only: real64


        integer :: i
        real(real64) :: result

        ! Integrate exp(-2x) from 0 to infinity. Expected result is 1/2.
        do i = 1, NNN1
            TA(i) = exp(-2*R(i)) * RP(i)
        end do
        call quad(result)

        print *, "RESULT from QUAD", result
        print *, "expected        ", 0.5_real64
        print *, "difference      ", abs(result - 0.5_real64)

        if(abs(result - 0.5_real64) > 1e-12_real64) then
            stop 1
        end if
    end subroutine runtest

    subroutine setup
        use grasp_parameters, only: nnn1, nnnp, radial_workspace_capacity
        use radial_grid_state, only: h => transformed_grid_step, hp => asymptotic_grid_step, n => radial_point_count, &
           & nnn1 => radial_workspace_capacity, r => radial_grid, rnt => radial_grid_scale, rp => radial_grid_derivative, &
           & rpor => radial_derivative_over_radius
        use radial_quadrature_workspace, only: nnn1 => radial_workspace_capacity, mtp, ta
        use radial_stencil_weights, only: setqic
        use radial_grid_setup, only: radgrd
        use, intrinsic :: iso_fortran_env, only: real64

        H = 5.0e-2_real64
        RNT = 2.0e-6_real64
        HP = 0.0_real64
        N = NNNP
        MTP = NNNP

        call setqic
        call radgrd
    end subroutine setup

end program core_quadrature
