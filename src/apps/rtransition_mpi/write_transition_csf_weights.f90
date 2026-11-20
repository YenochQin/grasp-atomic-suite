module write_transition_csf_weights_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_transition_csf_weights
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_transition_csf_weights &
        & (coulomb_or_magnetic_amplitude,babushkin_amplitude,initial_level_index,final_level_index)
        ! Module dependencies
        !                                                                      *
        !   This routine calculates  the CSF Coulomb, Babuskin, and magnetic   *
        !   matrix elements for  a transition  between  levels  separated by   *
        !   energy OMEGA.                                                      *
        !                                                                      *
        !   Modified for different initial and final state orbitals            *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: orbital_index_radix
        use biorthogonal_orbital_sets, only: final_csf_count, initial_csf_count
        use debug_output_settings, only: cutoff, rmcdhf_debug_flags
        use eigenstate_expansions, only: final_csf_mixing_coefficients, initial_csf_mixing_coefficients
        use angular_symmetry_labels, only: right_justified_j_labels, parity_labels
        use orbital_csf_basis, only: orbital_l, orbital_principal_n, orbital_symmetry_labels
        use transition_matrix_workspace, only: transition_multipole_rank, transition_multipole_type, &
            transition_angular_coefficients, transition_bra_csf_indices, transition_ket_csf_indices, &
            transition_integral_count, transition_integral_coefficient_offsets, packed_transition_integral_labels
        use selected_level_indices, only: final_level_indices, initial_level_indices
        use level_symmetry_metadata, only: final_level_parity_signs, initial_level_parity_signs, final_level_two_j_plus_one, &
            initial_level_two_j_plus_one

        ! Procedure interfaces
        use single_particle_transition_element_interface, only: single_particle_transition_element
        implicit none (type, external)

        ! Arguments
        real(real64), intent(out) :: coulomb_or_magnetic_amplitude, babushkin_amplitude
        integer, intent(in) :: initial_level_index, final_level_index

        ! Local constants
        integer, parameter :: nca=65536

        ! Local declarations
        integer :: ja, jb, m, nclr, ncup, ncount, nfi, ipar, jpar, &
            ia, ib, j
        real(real64) :: coulomb_single_particle_element, babushkin_single_particle_element, magnetic_single_particle_element, &
            & eps, mixing_coefficient_product, coeff, contr
        character(len=4) :: jlabi, jlabj, jpari, jparj
        !-----------------------------------------------------------------------
        !
        eps = 1.0e-10_real64
        !
        coulomb_or_magnetic_amplitude = 0.0_real64
        babushkin_amplitude = 0.0_real64
        ncount = 0
        nfi = 0
        !
        !Rasa -- start
        !
        !   J/pi labels for levels
        !
        jlabi = right_justified_j_labels(initial_level_two_j_plus_one(initial_level_index))
        jlabj = right_justified_j_labels(final_level_two_j_plus_one(final_level_index))
        ipar = (initial_level_parity_signs(initial_level_index)+3)/2
        jpar = (final_level_parity_signs(final_level_index)+3)/2
        jpari = parity_labels(ipar)
        jparj = parity_labels(jpar)

        if (rmcdhf_debug_flags(18)) then
            write (*,'("Upper level ",i3,a4,a4,10x," Lower level",i3,a4,a4)') &
                final_level_indices(final_level_index), jlabj, jparj, initial_level_indices(initial_level_index),    &
                jlabi, jpari
            write (*,'(4x,"i",3x,"Coeff.(i)",5x,"f",2x,"Coeff.(f)",2x,"Orb(i)",1x,"Orb(f)",2x,"MCT coeff.",4x,&
                &"Radial Factor",5x,"Contribution")')
        end if
        ! Rasa -- end
        do m = 1,transition_integral_count
            ja = packed_transition_integral_labels(m)
            jb = mod (ja,orbital_index_radix)
            ja = ja/orbital_index_radix
            nclr = transition_integral_coefficient_offsets(m)+1
            ncup = transition_integral_coefficient_offsets(m+1)
            if (nclr > ncup) cycle
            if (mod(orbital_l(ja)+orbital_l(jb)+transition_multipole_type+transition_multipole_rank,2) /= 0) cycle
            call single_particle_transition_element &
                & (ja,jb,coulomb_single_particle_element,babushkin_single_particle_element,magnetic_single_particle_element)
            !
            do j = nclr,ncup
                !
                ia = transition_bra_csf_indices(j)
                ib = transition_ket_csf_indices(j)
                !         IB = MOD (IA,NCA)
                !         IA = IA/NCA
                !
                !  Check on consistency
                !
                if (ia > initial_csf_count .or. ib > final_csf_count) cycle
                !
                !  Observ that IA and IB refer to the merged list
                !  whereas EVECII and EVECFF refers to the initial and
                !  final state lists. Therefore we have IB -> IB-NCFII
                !
                mixing_coefficient_product &
                    & =initial_csf_mixing_coefficients(ia+(initial_level_index-1)*initial_csf_count)*                          &
                    final_csf_mixing_coefficients(ib+(final_level_index-1)*final_csf_count)
                !     :           EVECFF(IB-NCFII+(LEV2-1)*NCFFF)
                coeff = transition_angular_coefficients(j)
                if (abs(mixing_coefficient_product) > eps) then
                    if (transition_multipole_type == 0) then
                        coulomb_or_magnetic_amplitude = &
                            & coulomb_or_magnetic_amplitude+coulomb_single_particle_element*coeff*mixing_coefficient_product
                        babushkin_amplitude = babushkin_amplitude+babushkin_single_particle_element*coeff*mixing_coefficient_product
                        !Rasa -- start
                        if (rmcdhf_debug_flags(18) .and. &
                            & (dabs(coulomb_single_particle_element*coeff*mixing_coefficient_product)          &
                            .gt. cutoff)) then
                            write (*,'(i5,2x,d12.5,i5,2x,d12.5,2(1x,1i2,1a2),1p,d15.6,1x,d15.6,1x,1a,1x,d15.6)') &
                                ia, initial_csf_mixing_coefficients(ia+(initial_level_index-1)*initial_csf_count),       &
                                ib, final_csf_mixing_coefficients(ib+(final_level_index-1)*final_csf_count),                 &
                                orbital_principal_n(ja),orbital_symmetry_labels(ja),orbital_principal_n(jb), &
                                & orbital_symmetry_labels(jb),coeff,coulomb_single_particle_element,"C",  &
                                coulomb_single_particle_element*coeff*mixing_coefficient_product
                            write(*,'(62x,1p,d15.6,1x,1a,1x,d15.6)') &
                                & babushkin_single_particle_element,"B", &
                                babushkin_single_particle_element*coeff*mixing_coefficient_product
                        end if
                    !Rasa -- end
                    else
                        contr=magnetic_single_particle_element*coeff*mixing_coefficient_product
                        coulomb_or_magnetic_amplitude = &
                            & coulomb_or_magnetic_amplitude+magnetic_single_particle_element*coeff*mixing_coefficient_product
                        !Rasa -- start
                        if (rmcdhf_debug_flags(18) .and. (abs(contr) > cutoff)) then
                            write (*,'(i5,2x,d12.5,i5,2x,d12.5,2(1x,1i2,1a2),1p,d15.6,1x,d15.6,1x,1a,1x,d15.6)') &
                                ia, initial_csf_mixing_coefficients(ia+(initial_level_index-1)*initial_csf_count),       &
                                ib, final_csf_mixing_coefficients(ib-(final_level_index-1)*final_csf_count),                 &
                                orbital_principal_n(ja),orbital_symmetry_labels(ja),orbital_principal_n(jb), &
                                & orbital_symmetry_labels(jb),coeff,magnetic_single_particle_element," ",   &
                                contr
                        end if
                        !Rasa -- end
                    end if
                end if
            end do
        end do
        !
        return
        !
        write (*,'("Upper level ",i3,a4,a4,10x," Lower level",i3,a4,a4)') &
            final_level_indices(final_level_index), jlabj, jparj, initial_level_index, jlabi, jpari
        !
    end subroutine write_transition_csf_weights
end module write_transition_csf_weights_interface
