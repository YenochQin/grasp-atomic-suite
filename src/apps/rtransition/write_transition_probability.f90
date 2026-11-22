! Report radiative rates, oscillator strengths, and line strengths for one transition.
module write_transition_probability_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_transition_probability
contains

    ! This  routine  prints the basic oscillator strength  information
    ! for transitions between level I and level J.
    !
    ! Last revision: 28 Dec 1992
    subroutine write_transition_probability(coulomb_or_magnetic_amplitude, babushkin_amplitude, i, j, &
        & transition_angular_frequency, factor, lines, same_initial_and_final_state)
        ! Module dependencies
        use calculation_file_units, only: summary_output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: einstein_a_au_to_si, einstein_b_au_to_si, pi, vacuum_speed_of_light_au
        use transition_output_cutoff, only: transition_probability_cutoff
        use angular_symmetry_labels, only: format_j_label, parity_labels
        use transition_matrix_workspace, only: transition_multipole_rank, transition_multipole_type, &
            total_coulomb_probabilities, total_babushkin_probabilities, transition_option_flags
        use level_symmetry_metadata, only: final_level_parity_signs, initial_level_parity_signs, final_level_two_j_plus_one, &
            initial_level_two_j_plus_one
        use selected_level_indices, only: final_level_indices, initial_level_indices
        implicit none (type, external)

        ! Arguments
        logical, intent(in) :: same_initial_and_final_state
        integer :: i
        integer, intent(in) :: j
        integer, intent(inout) :: lines
        real(real64), intent(in) :: coulomb_or_magnetic_amplitude, babushkin_amplitude
        real(real64), intent(in) :: transition_angular_frequency, factor

        ! Local declarations
        integer :: ipar, jpar, isigna, isignb
        real(real64) :: dblfac(10)
        real(real64) :: stfac, dll1, omc, faau, fosc, fbau, eng, ggfactor, &
            acsq, absq, ac, ab, bc, bb, oscc, oscb, sa, sb, ams, am, bm, oscm
        character(:), allocatable :: jlabi, jlabj
        character(len=4) :: jpari, jparj
        character(len=2) :: f1,f2

        ! Initial values
        !
        data dblfac/ 3.0000000000_real64, 1.5000000000e1_real64, 1.0500000000e2_real64, &
            9.4500000000e2_real64, 1.0395000000e4_real64, 1.3513500000e5_real64, 2.0270250000e6_real64, &
            3.4459425000e7_real64, 6.5472907500e8_real64, 1.3749310575e10_real64/
        ! Evaluate statistical factors and constants
        stfac = initial_level_two_j_plus_one(i)
        dll1 = dble(transition_multipole_rank + transition_multipole_rank + 1)
        stfac = stfac/dll1
        omc = transition_angular_frequency/vacuum_speed_of_light_au
        ! OMC = OMEGA/C
        faau = 2.0_real64*omc*stfac/dble(final_level_two_j_plus_one(j))
        fosc = vacuum_speed_of_light_au*stfac/omc
        ! FOSC = C*STFAC/OMC
        fbau = pi*fosc/transition_angular_frequency
        eng = transition_angular_frequency*factor
        if (transition_option_flags(1)) eng = factor/transition_angular_frequency
        ! J/pi labels for levels
        jlabi = format_j_label(initial_level_two_j_plus_one(i))
        jlabj = format_j_label(final_level_two_j_plus_one(j))
        ipar = (initial_level_parity_signs(i) + 3)/2
        jpar = (final_level_parity_signs(j) + 3)/2
        jpari = parity_labels(ipar)
        jparj = parity_labels(jpar)
        ggfactor = vacuum_speed_of_light_au**(2*transition_multipole_rank - 2)*dble(transition_multipole_rank)* &
            dblfac(transition_multipole_rank)**2/((2.0_real64*dble(transition_multipole_rank) + &
            1.0_real64)*(dble(transition_multipole_rank) + &
                & 1.0_real64)*abs(transition_angular_frequency)**(2*transition_multipole_rank - 1))
        ! end
        !
        !    Calculate Einstein A and B coefficients and oscillator strengths
        if (transition_multipole_type == 0) then
            ! Electric multipoles
            !
            ! In atomic units
            if (coulomb_or_magnetic_amplitude < 0) then
                isigna = -1
            else
                isigna = 1
            end if
            if (babushkin_amplitude < 0) then
                isignb = -1
            else
                isignb = 1
            end if

            acsq = coulomb_or_magnetic_amplitude**2
            absq = babushkin_amplitude**2
            ac = acsq*faau
            ab = absq*faau
            bc = acsq*fbau
            bb = absq*fbau
            oscc = acsq*fosc
            oscb = absq*fosc
            sa = oscc*ggfactor
            sb = oscb*ggfactor
            ! Convert the transition rates and coefficients to SI when requested.
            if (.not.transition_option_flags(7)) then
                ac = ac*einstein_a_au_to_si
                ab = ab*einstein_a_au_to_si
                bc = bc*einstein_b_au_to_si
                bb = bb*einstein_b_au_to_si
            end if
            ! Accumulate total of A coefficients
            total_coulomb_probabilities(j) = total_coulomb_probabilities(j) + ac
            total_babushkin_probabilities(j) = total_babushkin_probabilities(j) + ab
            ! Print information if both AC and AB are greater than CUTOFF
            if (abs(ac)>=transition_probability_cutoff .and. abs(ab)>=transition_probability_cutoff) then
                if (-1.0e9_real64 < eng .and. eng < -1.0e-2_real64) then
                    if (same_initial_and_final_state) then
                        f1 = "f "
                        f2 = "f "
                    else
                        f1 = "f2"
                        f2 = "f1"
                    end if
                    write (summary_output_unit, '(1x,a2,i3,1x,2a,a2,i3,1x,2a,0p,f13.2," C",1p,3d13.5)') &
                        f1,final_level_indices(j),jlabj,jparj,f2,initial_level_indices(i),jlabi,jpari, &
                        -eng, -ac, -oscc, -sa !, ASFA
                    write (summary_output_unit, '(42x," B",1p,3d13.5)') -ab, -oscb, -sb !, ASFB
                else if ( 1.0e9_real64 > eng .and. eng > 1.0e-2_real64) then
                    if (same_initial_and_final_state) then
                        f1 = "f "
                        f2 = "f "
                    else
                        f1 = "f1"
                        f2 = "f2"
                    end if
                    write (summary_output_unit,'(1x,a2,i3,1x,2a,a2,i3,1x,2a,0p,f13.2," C",1p,3d13.5)') &
                        f1,initial_level_indices(i),jlabi,jpari,f2,final_level_indices(j), &
                        jlabj,jparj,eng,ac*final_level_two_j_plus_one(j)/initial_level_two_j_plus_one(i),     &
                        oscc,sa
                    write (summary_output_unit,'(42x," B",1p,3d13.5)') ab* &
                        final_level_two_j_plus_one(j)/initial_level_two_j_plus_one(i),oscb,sb
                end if
                lines = lines + 3
            end if
        else
            ! Magnetic multipoles
            !
            ! In atomic units
            if (coulomb_or_magnetic_amplitude < 0) then
                isigna = -1
            else
                isigna = 1
            end if
            ams = coulomb_or_magnetic_amplitude**2
            am = ams*faau
            bm = ams*fbau
            oscm = ams*fosc
            sa = oscm*ggfactor*vacuum_speed_of_light_au*vacuum_speed_of_light_au*4
            ! end
            !
            ! Convert the transition rates and coefficients to SI when requested.
            if (.not.transition_option_flags(7)) then
                am = am*einstein_a_au_to_si
                bm = am*einstein_b_au_to_si
            end if
            ! Accumulate total of A coefficients
            total_coulomb_probabilities(j) = total_coulomb_probabilities(j) + am
            total_babushkin_probabilities(j) = total_babushkin_probabilities(j) + am
            ! Print information if AM is greater than CUTOFF
            if (abs(am) >= transition_probability_cutoff) then
                if (-1.0e9_real64 < eng .and. eng < -1.0e-2_real64) then
                    if (same_initial_and_final_state) then
                        f1 = "f "
                        f2 = "f "
                    else
                        f1 = "f2"
                        f2 = "f1"
                    end if
                    write (summary_output_unit,'(1x,a2,i3,1x,2a,a2,i3,1x,2a,0p,f13.2," M",1p,3d13.5)') &
                        f1,final_level_indices(j),jlabj,jparj,f2,initial_level_indices(i),&
                        jlabi,jpari,-eng,-am,-oscm,-sa
                ! :               ,ASFA (M-value)
                else if ( 1.0e9_real64 > eng .and. eng > 1.0e-2_real64) then
                    if (same_initial_and_final_state) then
                        f1 = "f "
                        f2 = "f "
                    else
                        f1 = "f1"
                        f2 = "f2"
                    end if
                    write (summary_output_unit,'(1x,a2,i3,1x,2a,a2,i3,1x,2a,0p,f13.2," M",1p,3d13.5)') &
                        f1,initial_level_indices(i),jlabi,jpari,f2,final_level_indices(j),&
                        jlabj,jparj,eng,am*final_level_two_j_plus_one(j)/initial_level_two_j_plus_one(i),    &
                        oscm,sa
                    ! :              ,ASFA (M-value)
                end if
                lines = lines+2
            end if
        end if
        return
        ! format for highly charged ions         F11.2 -> F13.2
        ! format for highly charged ions         F11.2 -> F13.2
        return
    end subroutine write_transition_probability
end module write_transition_probability_interface
