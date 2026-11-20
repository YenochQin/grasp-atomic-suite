module transverse_radial_integral_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: transverse_radial_integral
contains

    !***********************************************************************
    !                                                                      *
    real(real64) function transverse_radial_integral (recoupling_coefficient, orbital_a, orbital_c, rbd, orbital_b, &
        & orbital_d, multipole_rank, iw)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        !                                                                      *
        !   This routine evaluates transverse interaction integrals:           *
        !                                                                      *
        !                          (k)                                         *
        !                         S   (a,c;b,d;w)                              *
        !                                                                      *
        !   where w = wac if IW = 1, and w= wbd if IW = 2.                     *
        !                                                                      *
        !   Call(s) to: [LIB92]: QUAD.                                         *
        !               [RCI92]: ZKF.                                          *
        !                                           Last update: 06 Nov 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_radial_points
        use bessel_function_workspace, only: breit_spherical_bessel_values, breit_spherical_neumann_values, &
            max_radial_points, breit_bessel_frequencies
        use debug_output_settings, only: rmcdhf_debug_flags
        use radial_grid_state, only: radial_grid, radial_grid_derivative, radial_derivative_over_radius
        use orbital_csf_basis, only: orbital_principal_n, orbital_symmetry_labels
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a, radial_workspace_b
        use radial_wavefunctions, only: orbital_last_grid_index

        ! Procedure interfaces
        use evaluate_breit_hartree_z_interface, only: evaluate_breit_hartree_z

        ! Module dependencies
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_a, orbital_b, orbital_c, orbital_d
        integer, intent(in) :: multipole_rank, iw
        real(real64), intent(in) :: recoupling_coefficient(max_radial_points)
        real(real64), intent(in) :: rbd(max_radial_points)

        ! Local declarations
        integer :: mxrac, mxrbd
        real(real64) :: tkeep(max_radial_points)
        real(real64) :: epsi, w, wk, valu
        !-----------------------------------------------
        !
        !
        epsi = 1.0e-10_real64
        !
        w = breit_bessel_frequencies(iw)
        wk = dble(multipole_rank + multipole_rank + 1)
        mxrac = min(orbital_last_grid_index(orbital_a),orbital_last_grid_index(orbital_c))
        mxrbd = min(orbital_last_grid_index(orbital_b),orbital_last_grid_index(orbital_d))
        !
        !           (k-1)
        !  Compute Z     (rho  ; s)
        !                    ac
        !
        radial_workspace_a(:mxrac) = recoupling_coefficient(:mxrac)
        integration_last_index = mxrac
        call evaluate_breit_hartree_z (multipole_rank - 1, orbital_a, orbital_c)
        !
        if (dabs(w) < epsi) then
            !
            !   W = 0 case
            !
            tkeep(:integration_last_index) = radial_workspace_b(:integration_last_index)
            call evaluate_breit_hartree_z (multipole_rank + 1, orbital_a, orbital_c)
            integration_last_index = min(integration_last_index,mxrbd)
            radial_workspace_a(1) = 0.0_real64
            radial_workspace_a(2:integration_last_index) = &
                rbd(2:integration_last_index)*radial_derivative_over_radius(2:integration_last_index)* &
                    (radial_workspace_b(2:integration_last_index)-tkeep(2:integration_last_index))
            call integrate_radial_grid (valu)
            transverse_radial_integral = wk*valu*0.5_real64
        !
        else
            !
            !   Finite w: see I P Grant and B J McKenzie, J Phys B: At Mol Phys,
            !   13 (1980) 2671-2681
            !
            tkeep(:integration_last_index) = radial_workspace_b(:integration_last_index)
            radial_workspace_a(:integration_last_index) = -radial_workspace_a(:integration_last_index)* &
                breit_spherical_bessel_values(1,iw, &
                :integration_last_index)
            call evaluate_breit_hartree_z (multipole_rank - 1, orbital_a, orbital_c)
            !
            integration_last_index = min(integration_last_index,mxrbd)
            radial_workspace_a(1) = 0.0_real64
            radial_workspace_a(2:integration_last_index) = ((1.0_real64 + breit_spherical_neumann_values(2,iw, &
                2:integration_last_index))*radial_workspace_b(2:integration_last_index)-tkeep(2:integration_last_index)* &
                    breit_spherical_neumann_values&
                (2,iw, &
                    2:integration_last_index))*rbd(2:integration_last_index)/radial_grid(2:integration_last_index)**2* &
                        radial_derivative_over_radius(2:integration_last_index)
            call integrate_radial_grid (valu)
            transverse_radial_integral = (wk/w)**2*valu
            radial_workspace_a(:mxrbd) = rbd(:mxrbd)*(1.0_real64 + breit_spherical_bessel_values(2,iw,:mxrbd))
            integration_last_index = mxrbd
            call evaluate_breit_hartree_z (multipole_rank + 1, orbital_b, orbital_d)
            integration_last_index = min(integration_last_index,mxrac)
            radial_workspace_a(1) = 0.0_real64
            radial_workspace_a(2:integration_last_index) = recoupling_coefficient(2:integration_last_index)*(1.0_real64 + &
                breit_spherical_neumann_values(1,iw, &
                2:integration_last_index))*radial_workspace_b(2:integration_last_index)*radial_grid(2:integration_last_index)&
                *radial_grid_derivative(2:integration_last_index)
            call integrate_radial_grid (valu)
            transverse_radial_integral = transverse_radial_integral - valu*w*w/dble((2*multipole_rank + 3)*(2*multipole_rank - 1))
            !
        end if
        !
        if (rmcdhf_debug_flags(11)) then
            call ensure_debug_output()
            write (angular_debug_unit, '("  (",1i2,")"/,"S     (",1i2,1a2,",",1i2,1a2,"|",1i2,1a2,",",1i2,1a2,";",1i2,&
                &") = ",1p,d19.12)') multipole_rank, orbital_principal_n(orbital_a), orbital_symmetry_labels(orbital_a), &
                & orbital_principal_n(orbital_c), &
                & orbital_symmetry_labels(orbital_c), orbital_principal_n(orbital_b)&
                , orbital_symmetry_labels(orbital_b), orbital_principal_n(orbital_d), orbital_symmetry_labels(orbital_d), iw, valu
        end if
        !
        return
        !
        return
        !
    end function transverse_radial_integral
end module transverse_radial_integral_interface
