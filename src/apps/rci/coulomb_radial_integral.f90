module coulomb_radial_integral_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: coulomb_radial_integral
contains

    !***********************************************************************
    !                                                                      *
    real(real64) function coulomb_radial_integral (recoupling_coefficient, orbital_a, orbital_c, rbd, orbital_b, orbital_d, &
        & multipole_rank, iw)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        !                                                                      *
        !   This routine evaluates the transverse interaction integrals.  If   *
        !   IW = 0, it calulates the U(r1,r2) integral; otherwise, it calcu-   *
        !   lates  R bar (k; a c | b d ; w)  with  w = wac  if  IW = 1,  and   *
        !   w = wbd  if IW = 2.                                                *
        !                                                                      *
        !   Call(s) to: [LIB92]: QUAD.                                         *
        !               [RCI92]: ZKF.                                          *
        !                                           Last update: 15 Oct 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_radial_points
        use bessel_function_workspace, only: breit_spherical_bessel_values, breit_spherical_neumann_values, max_radial_points
        use debug_output_settings, only: rmcdhf_debug_flags
        use radial_grid_state, only: radial_derivative_over_radius
        use orbital_csf_basis, only: orbital_principal_n, orbital_symmetry_labels
        use radial_quadrature_workspace,  only: integration_last_index, radial_workspace_a, radial_workspace_b
        use radial_wavefunctions, only: orbital_last_grid_index

        ! Procedure interfaces
        use evaluate_breit_hartree_z_interface, only: evaluate_breit_hartree_z

        ! Module dependencies
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_a, orbital_c
        integer, intent(in) :: orbital_b, orbital_d
        integer, intent(in) :: multipole_rank
        integer, intent(in) :: iw
        real(real64), intent(in) :: recoupling_coefficient(max_radial_points)
        real(real64), intent(in) :: rbd(max_radial_points)

        ! Local declarations
        integer :: mxrbd, mxrac
        real(real64) :: result
        !-----------------------------------------------
        !
        mxrbd = min (orbital_last_grid_index(orbital_b),orbital_last_grid_index(orbital_d))
        mxrac = min (orbital_last_grid_index(orbital_a),orbital_last_grid_index(orbital_c))
        !
        if (iw == 0) then
            !
            !   IW = 0
            !
            radial_workspace_a(:mxrac) = recoupling_coefficient(:mxrac)
            integration_last_index = mxrac
            call evaluate_breit_hartree_z (multipole_rank, orbital_a, orbital_c)
            integration_last_index = min(integration_last_index,mxrbd)
            radial_workspace_a(1) = 0.0_real64
            radial_workspace_a(2:integration_last_index) = &
                rbd(2:integration_last_index)*radial_workspace_b(2:integration_last_index)* &
                    radial_derivative_over_radius(2:integration_last_index)
        !
        else
            !
            !   IW = 1,2
            !
            radial_workspace_a(:mxrac) = recoupling_coefficient(:mxrac)*(1.0_real64 + breit_spherical_bessel_values(1,iw,:mxrac))
            integration_last_index = mxrac
            call evaluate_breit_hartree_z (multipole_rank, orbital_a, orbital_c)
            integration_last_index = min(integration_last_index,mxrbd)
            radial_workspace_a(1) = 0.0_real64
            radial_workspace_a(2:integration_last_index) = rbd(2:integration_last_index)*(1.0_real64 + &
                breit_spherical_neumann_values(1,iw, &
                2:integration_last_index))*radial_workspace_b(2:integration_last_index)*radial_derivative_over_radius(2:&
                integration_last_index)
            !
        end if
        !
        call integrate_radial_grid (result)
        coulomb_radial_integral = result
        !
        !   Debug printout if option set
        !
        if (rmcdhf_debug_flags(11)) then
            call ensure_debug_output()
            write (angular_debug_unit, '("_ (",1i2,")"/,"R     (",1i2,1a2,",",1i2,1a2,"|",1i2,1a2,",",1i2,1a2,";",1i2,&
                &") = ",1p,d19.12)') multipole_rank, orbital_principal_n(orbital_a), orbital_symmetry_labels(orbital_a), &
                & orbital_principal_n(orbital_c), &
                & orbital_symmetry_labels(orbital_c), orbital_principal_n(orbital_b)&
                , orbital_symmetry_labels(orbital_b), orbital_principal_n(orbital_d), orbital_symmetry_labels(orbital_d), iw, result
        end if
        !
        return
        !
        return
        !
    end function coulomb_radial_integral
end module coulomb_radial_integral_interface
