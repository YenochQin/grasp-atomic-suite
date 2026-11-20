module single_particle_transition_element_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: single_particle_transition_element
contains

    !***********************************************************************
    !                                                                      *
    subroutine single_particle_transition_element(first_orbital, second_orbital, coulomb_element, babushkin_element, &
        & magnetic_element)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        !                                                                      *
        !   This routine calculates the reduced matrix elements for pair I,J   *
        !   in either Coulomb/Babuskin gauge or for magnetic case.             *
        !                                                                      *
        !   These are defined in the  Brink and Satchler sense --- i.e. com-   *
        !   patible with Pyper et al.'s MCT paper but not with Grant's radi-   *
        !   ative transitions paper.                                           *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use bessel_function_workspace, only: transition_spherical_bessel_values, cross_component_product_difference, &
            same_component_product_sum
        use biorthogonal_orbital_sets, only: merged_to_final_orbital_index, merged_to_initial_orbital_index
        use debug_output_settings, only: rmcdhf_debug_flags
        use atomic_calculation_state, only: speed_of_light_au, vacuum_speed_of_light_au
        use radial_grid_state, only: radial_point_count, radial_grid_derivative
        use orbital_csf_basis, only: orbital_kappa, orbital_two_j
        use transition_matrix_workspace, only: transition_multipole_rank, transition_multipole_type
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a, radial_workspace_b
        use radial_wavefunctions, only: final_orbital_last_grid_index, initial_orbital_last_grid_index, &
            final_radial_large_components, final_radial_small_components, initial_radial_large_components, &
            initial_radial_small_components
        use relativistic_three_j, only: relativistic_three_j_value
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: first_orbital,second_orbital
        real(real64), intent(out) :: coulomb_element, babushkin_element, magnetic_element

        ! Local declarations
        integer :: nkji, nkjj, ipj, naki, nakj, ii, lp, lm
        real(real64) :: eps, gauge_correction, tji, tjj, fact, form, fl, flp, dfki, dfkj, &
            upper_order_coulomb_scale, lower_order_coulomb_scale, upper_order_cross_sum_factor, &
                & upper_order_cross_diff_factor, lower_order_cross_sum_factor, lower_order_cross_diff_factor, value, cjl, &
                & cip, tarom, esther, &
            tekel, peres, tamar, enoch, seth, shem, darius, cyrus, daniel, boaz, &
            esau, aaron, gimel
        !-----------------------------------------------

        !
        eps = 1.0e-10_real64
        !
        if (rmcdhf_debug_flags(12) .or. rmcdhf_debug_flags(13) .or. rmcdhf_debug_flags(14) .or. rmcdhf_debug_flags(15)) then
            call ensure_debug_output()
            write (angular_debug_unit, &
                '(/,1x,30("+"),"   SPME called for orbitals ",i4,"   and",i4,3x,30("+"))') first_orbital, second_orbital
        end if
        !
        magnetic_element = 0.0_real64
        coulomb_element = 0.0_real64
        gauge_correction = 0.0_real64
        babushkin_element = 0.0_real64
        !
        !  Evaluate factor multiplying Mbar(a,b)
        !
        nkji = orbital_two_j(first_orbital)
        nkjj = orbital_two_j(second_orbital)
        tji = dble(nkji)
        tjj = dble(nkjj)
        fact = 1.0_real64
        ipj = (nkji + 1)/2 + nkjj
        if (mod(ipj,2) /= 0) fact = -fact
        naki = orbital_kappa(first_orbital)
        nakj = orbital_kappa(second_orbital)
        form = fact*sqrt(tjj + 1.0_real64)*relativistic_three_j_value(nakj,transition_multipole_rank,naki)
        fl = dble(transition_multipole_rank)
        flp = fl + 1.0_real64
        dfki = dble(naki)
        dfkj = dble(nakj)
        !
        integration_last_index = min(initial_orbital_last_grid_index(merged_to_initial_orbital_index(first_orbital)), &
            final_orbital_last_grid_index(merged_to_final_orbital_index(second_orbital)))
        !
        if (transition_multipole_type == 0) then
            !
            if (abs(form) > eps) then
                !
                !  To pick the right initial and final state radial functions convert
                !  from orbital order of the merged list to the orbital orders of the
                !  initial and final state lists, respectively.
                !  N(I) -> NNII(I), N(J) -> NNFF(J)
                !
                do ii = 1, radial_point_count
                    radial_workspace_b(ii) = initial_radial_large_components(ii, &
                        merged_to_initial_orbital_index(first_orbital))*final_radial_small_components(ii, &
                        merged_to_final_orbital_index(second_orbital)) + initial_radial_small_components(ii, &
                        merged_to_initial_orbital_index(first_orbital))*&
                        final_radial_large_components(ii,merged_to_final_orbital_index(second_orbital))
                    cross_component_product_difference(ii) = initial_radial_large_components(ii, &
                        merged_to_initial_orbital_index(first_orbital))* &
                        final_radial_small_components(ii,merged_to_final_orbital_index(second_orbital)) - &
                        initial_radial_small_components(ii,merged_to_initial_orbital_index(first_orbital))*&
                        final_radial_large_components(ii,merged_to_final_orbital_index(second_orbital))
                    ! zou
                    radial_workspace_b(ii) = radial_workspace_b(ii)*speed_of_light_au/vacuum_speed_of_light_au
                    cross_component_product_difference(ii) = cross_component_product_difference(ii)* &
                        speed_of_light_au/vacuum_speed_of_light_au
                    same_component_product_sum(ii) = initial_radial_large_components(ii, &
                        merged_to_initial_orbital_index(first_orbital))* &
                        final_radial_large_components(ii,merged_to_final_orbital_index(second_orbital)) + &
                        initial_radial_small_components(ii,merged_to_initial_orbital_index(first_orbital))*&
                        final_radial_small_components(ii,merged_to_final_orbital_index(second_orbital))* &
                            (speed_of_light_au/vacuum_speed_of_light_au)**2
                    !           TD(II) = PFII(II,NNII(I))*PFFF(II,NNFF(J)) +
                    !    :               QFII(II,NNII(I))*QFFF(II,NNFF(J))
                    ! zou
                end do
                lp = transition_multipole_rank + 1
                lm = transition_multipole_rank - 1
                !
                !  Calculate Coulomb coefficients
                !
                upper_order_coulomb_scale = sqrt(fl/flp)
                lower_order_coulomb_scale = -sqrt(flp/fl)
                upper_order_cross_sum_factor = upper_order_coulomb_scale*(dfki - dfkj)
                upper_order_cross_diff_factor = upper_order_coulomb_scale*flp
                lower_order_cross_sum_factor = lower_order_coulomb_scale*(dfki - dfkj)
                lower_order_cross_diff_factor = -lower_order_coulomb_scale*fl
                !
                !  Tabulate Coulomb integrand
                !
                radial_workspace_a(1) = 0.0_real64
                radial_workspace_a(2:integration_last_index) = &
                    radial_grid_derivative(2:integration_last_index)*(transition_spherical_bessel_values(2:integration_last_index, &
                    3)*(upper_order_cross_sum_factor*radial_workspace_b(2:integration_last_index)+upper_order_cross_diff_factor* &
                        cross_component_product_difference(2:integration_last_index)&
                    )+transition_spherical_bessel_values(2:integration_last_index,1)*(lower_order_cross_sum_factor* &
                        radial_workspace_b(2:integration_last_index)+lower_order_cross_diff_factor* &
                        cross_component_product_difference(2:integration_last_index)))
                call integrate_radial_grid (value)
                coulomb_element = form*value
                !
                !  Calculate gauge dependent coefficients
                !
                cjl = -(fl + flp)
                cip = dfki - dfkj
                !
                !  Tabulate gauge dependent integrand
                !
                radial_workspace_a(1) = 0.0_real64
                radial_workspace_a(2:integration_last_index) = &
                    radial_grid_derivative(2:integration_last_index)* &
                        (transition_spherical_bessel_values(2:integration_last_index,2)*(cjl* &
                        same_component_product_sum(2:integration_last_index))+ &
                    transition_spherical_bessel_values(2:integration_last_index,3)*(cip&
                    *radial_workspace_b(2:integration_last_index)+flp* &
                        cross_component_product_difference(2:integration_last_index))+ &
                        transition_spherical_bessel_values(2:integration_last_index, &
                        1)*(cip*radial_workspace_b(2:integration_last_index)-fl* &
                            cross_component_product_difference(2:integration_last_index&
                    )))
                !
                !  Print gauge dependent integrand if requested
                !
                if (rmcdhf_debug_flags(13)) then
                    write (angular_debug_unit, '(/,1x,20x,"Local form of gauge dependent integral for"," orbitals",i3,&
                        &" and",i3,/,/,100(1x,7(i3,2x,1p,d11.3,2x),/))') first_orbital, second_orbital, (ii, &
                            radial_workspace_a(ii),ii=1, &
                        & radial_point_count)
                end if
                call integrate_radial_grid (value)
                gauge_correction = form*value
                babushkin_element = coulomb_element + sqrt(flp/fl)*gauge_correction
                !
                !  Print Coulomb and gauge dependent integrals if requested
                !
                if (rmcdhf_debug_flags(12)) then
                    write (angular_debug_unit, '(/,1x,"Orbital pair (",i2,",",i2,") integrals:","Coul. gauge ",1p,d13.3,3x,&
                        &"Gauge contribution",d13.3,3x,"Bab.gauge",d13.3)') first_orbital, second_orbital, coulomb_element, &
                            & gauge_correction, babushkin_element
                end if
                !
                !  Calculate the contibutions from various terms if requested
                !
                if (rmcdhf_debug_flags(14) .or. rmcdhf_debug_flags(15)) then
                    !
                    radial_workspace_a(1) = 0.0_real64
                    radial_workspace_a(2:integration_last_index) = same_component_product_sum(2:integration_last_index)* &
                        transition_spherical_bessel_values(2:integration_last_index, &
                        2)*radial_grid_derivative(2:integration_last_index)
                    if (rmcdhf_debug_flags(15)) write (angular_debug_unit, '(/,100(1x,7(i3,2x,1p,d11.3,2x),/),/)') (ii, &
                        radial_workspace_a(ii),ii=1, &
                        & radial_point_count)
                    call integrate_radial_grid (tarom)
                    if (rmcdhf_debug_flags(14)) write (angular_debug_unit, '(/," Integral J(L)  =  ",1p,d28.12)') tarom
                    esther = -cjl
                    if (rmcdhf_debug_flags(14)) write (angular_debug_unit, '(" Factor multiplying this ",1p,d20.12)') esther
                    tekel = esther*tarom
                    if (rmcdhf_debug_flags(14)) write (angular_debug_unit, '(" Resultant contribution  ",1p,d20.12,/)') tekel
                    !
                    radial_workspace_a(1) = 0.0_real64
                    radial_workspace_a(2:integration_last_index) = &
                        -radial_workspace_b(2:integration_last_index)*transition_spherical_bessel_values(2:integration_last_index, &
                        3)*radial_grid_derivative(2:integration_last_index)
                    if (rmcdhf_debug_flags(15)) write (angular_debug_unit, '(/,100(1x,7(i3,2x,1p,d11.3,2x),/),/)') (ii, &
                        radial_workspace_a(ii),ii=1, &
                        & radial_point_count)
                    call integrate_radial_grid (peres)
                    if (rmcdhf_debug_flags(14)) then
                        write (angular_debug_unit, '(" I + (L+1) integral  = ",1p,d24.12)') peres
                        write (angular_debug_unit, '(" Factor multiplying this ",1p,d20.12)') cip
                    end if
                    tamar = cip*peres
                    if (rmcdhf_debug_flags(14)) write (angular_debug_unit, '(" Resultant contribution  ",1p,d20.12,/)') tamar
                    !
                    radial_workspace_a(1) = 0.0_real64
                    radial_workspace_a(2:integration_last_index) = &
                        -radial_workspace_b(2:integration_last_index)*transition_spherical_bessel_values(2:integration_last_index, &
                        1)*radial_grid_derivative(2:integration_last_index)
                    if (rmcdhf_debug_flags(15)) write (angular_debug_unit, '(/,100(1x,7(i3,2x,1p,d11.3,2x),/),/)') (ii, &
                        radial_workspace_a(ii),ii=1, &
                        & radial_point_count)
                    call integrate_radial_grid (enoch)
                    if (rmcdhf_debug_flags(14)) then
                        write (angular_debug_unit, '(" I + (L-1) integral  = ",1p,d24.12)') enoch
                        write (angular_debug_unit, '(" Factor multiplying this ",1p,d20.12)') cip
                    end if
                    seth = cip*enoch
                    if (rmcdhf_debug_flags(14)) write (angular_debug_unit, '(" Resultant contribution  ",1p,d20.12,/)') seth
                    !
                    radial_workspace_a(1) = 0.0_real64
                    radial_workspace_a(2:integration_last_index) = transition_spherical_bessel_values(2:integration_last_index, &
                        1)*cross_component_product_difference(2:integration_last_index)* &
                            radial_grid_derivative(2:integration_last_index)
                    if (rmcdhf_debug_flags(15)) write (angular_debug_unit, '(/,100(1x,7(i3,2x,1p,d11.3,2x),/),/)') (ii, &
                        radial_workspace_a(ii),ii=1, &
                        & radial_point_count)
                    call integrate_radial_grid (shem)
                    if (rmcdhf_debug_flags(14)) then
                        write (angular_debug_unit, '(" I - (L-1) integral  = ",1p,d24.12)') shem
                        write (angular_debug_unit, '(" Factor multiplying this ",1p,d20.12)') fl
                    end if
                    darius = fl*shem
                    if (rmcdhf_debug_flags(14)) write (angular_debug_unit, '(" Resultant contribution  ",1p,d20.12,/)') darius
                    !
                    radial_workspace_a(1) = 0.0_real64
                    radial_workspace_a(2:integration_last_index) = -transition_spherical_bessel_values(2:integration_last_index, &
                        3)*cross_component_product_difference(2:integration_last_index)* &
                            radial_grid_derivative(2:integration_last_index)
                    if (rmcdhf_debug_flags(15)) write (angular_debug_unit, '(/,100(1x,7(i3,2x,1p,d11.3,2x),/),/)') (ii, &
                        radial_workspace_a(ii),ii=1, &
                        & radial_point_count)
                    call integrate_radial_grid (cyrus)
                    if (rmcdhf_debug_flags(14)) write (angular_debug_unit, '(" I - (L+1) integral  = ",1p,d24.12)') cyrus
                    daniel = flp
                    if (rmcdhf_debug_flags(14)) write (angular_debug_unit, '(" Factor multiplying this ",1p,d20.12)') daniel
                    boaz = daniel*cyrus
                    if (rmcdhf_debug_flags(14)) write (angular_debug_unit, '(" Resultant contribution  ",1p,d20.12,/)') boaz
                    !
                    esau = tamar + seth + darius + boaz
                    if (rmcdhf_debug_flags(14)) then
                        write (angular_debug_unit, '(" These last four add up to ",1p,d20.12,&
                            &"    compared to the first which is ",d20.12)') esau, tekel
                    end if
                    aaron = esau + tekel
                    gimel = -value
                    if (rmcdhf_debug_flags(14)) then
                        write (angular_debug_unit, '(" These all add up to ",1p,d20.12," which should equal ",d20.12)') aaron, gimel
                        write (angular_debug_unit, '(/,18x,"Integral",11x,"Coefficient",8x,"Contribution")')
                        write (angular_debug_unit, '(" I + (L+1)",1p,3d20.10)') peres, cip, tamar
                        write (angular_debug_unit, '(" I + (L-1)",1p,3d20.10)') enoch, cip, seth
                        write (angular_debug_unit, '(" I - (L-1)",1p,3d20.10)') shem, fl, darius
                        write (angular_debug_unit, '(" I - (L+1)",1p,3d20.10)') cyrus, daniel, boaz
                        write (angular_debug_unit, '(54x,16("-"),/,50x,1p,d20.10)') esau
                        write (angular_debug_unit, '(" J (L)",4x,1p,3d20.10)') tarom, esther, tekel
                        write (angular_debug_unit, '(54x,16("-"),/," Correct figure",1p,d20.10,15x,d20.10)') gimel, aaron
                    end if
                    !
                end if
                !
            end if
            !
            if (rmcdhf_debug_flags(12)) then
                write (angular_debug_unit, '(/,1x,"Orbital pair (",i2,",",i2,") integrals:","Coul. gauge ",1p,d13.3,3x,&
                    &"Gauge contribution",d13.3,3x,"Bab.gauge",d13.3)') first_orbital, second_orbital, coulomb_element, &
                        & gauge_correction, babushkin_element
            end if
        !
        else
            !
            !   Tabulate magnetic integrand
            !
            if (abs(form) > eps) then
                !
                radial_workspace_a(1) = 0.0_real64
                do ii = 2, integration_last_index
                    !zou
                    radial_workspace_a(ii) = (initial_radial_large_components(ii, &
                        merged_to_initial_orbital_index(first_orbital))*final_radial_small_components(ii, &
                        merged_to_final_orbital_index(second_orbital))+initial_radial_small_components(ii, &
                        merged_to_initial_orbital_index(first_orbital))*&
                        final_radial_large_components(ii,merged_to_final_orbital_index(second_orbital)))* &
                            transition_spherical_bessel_values(ii,2)* &
                            radial_grid_derivative(ii) &
                        & *speed_of_light_au/vacuum_speed_of_light_au
                    !           TA(II) =  (PFII(II,NNII(I))*QFFF(II,NNFF(J)) +
                    !    :      QFII(II,NNII(I))*PFFF(II,NNFF(J)))*BJ(II,2)*RP(II)
                    !zou
                end do
                call integrate_radial_grid (value)
                if (rmcdhf_debug_flags(14)) write (angular_debug_unit, '(" I + (L)  integral",1p,d20.10)') value
                magnetic_element = -value*(fl + flp)*(dfki + dfkj)*form/sqrt(fl*flp)
                !
            end if
            !
            if (rmcdhf_debug_flags(12)) write (angular_debug_unit, &
                '(/,1x,"Magnetic single particle matrix element ",1p,d13.3)') magnetic_element
            !
        end if
        !
        return
        !
        return
        !
    end subroutine single_particle_transition_element
end module single_particle_transition_element_interface
