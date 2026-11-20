!***********************************************************************
!                                                                      *
real(real64)  function vacuum_polarization_integral (first_orbital,second_orbital)
    ! Module dependencies
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    !                                                                      *
    !   Computes nuclear vacuum polarization integrals.                    *
    !                                                                      *
    !   Call(s) to: [LIB92]: QUAD.                                         *
    !                                           Last update: 09 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    !...Translated by Charlotte Froese Fischer
    !                       Gediminas Gaigalas  10/05/17
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: orbital_index_radix
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use nuclear_charge_distribution, only: max_radial_points, zdist
    use orbital_csf_basis, only: max_orbitals, orbital_principal_n, orbital_symmetry_labels
    use radial_quadrature_workspace, only: mtp, ta, tb
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    use radial_quadrature, only: quad
    implicit none (type, external)

    ! Arguments
    integer , intent(in) :: first_orbital, second_orbital

    ! Local constants
    integer, parameter :: orbital_label_radix = orbital_index_radix

    ! Local declarations
    integer :: k
    !------------------------------------------------
    !
    mtp = min (orbital_last_grid_index(first_orbital),orbital_last_grid_index(second_orbital))
    !TA(1) = 0.0D 00  !! Why would this be wrong?!!
    ta(1)  = 0
    do k = 2,mtp
        ta(k) = (radial_large_components(k,first_orbital)*radial_large_components(k,second_orbital)+radial_small_components(k, &
            & first_orbital)* &
            & radial_small_components(k,second_orbital))                  &
            *zdist(k)
    end do
    call quad (vacuum_polarization_integral)
    !
    if (ldbpr(9)) then
        call ensure_debug_output()
        write (angular_debug_unit,'(/"VPINTF: V (",1i2,1a2,",",1i2,1a2,") = ",1pd19.12)') &
            orbital_principal_n(first_orbital),orbital_symmetry_labels(first_orbital),orbital_principal_n(second_orbital), &
            & orbital_symmetry_labels(second_orbital),vacuum_polarization_integral
    end if
    !
    return
    !
    !
end function vacuum_polarization_integral
