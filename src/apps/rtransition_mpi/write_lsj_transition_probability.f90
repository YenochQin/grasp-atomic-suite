module write_lsj_transition_probability_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_lsj_transition_probability
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_lsj_transition_probability &
        & (inum_ii,inum_ff,coulomb_or_magnetic_amplitude,babushkin_amplitude,i,j,transition_angular_frequency,factor)
        ! Module dependencies
        use calculation_file_units, only: lsj_transition_output_unit
        use, intrinsic :: iso_fortran_env, only: real32
        !                                                                      *
        !   This  routine  prints the basic oscillator strength  information   *
        !   for transitions between level I and level J in LSJ coupling.       *
        !                                                                      *
        !   Written by G. Gaigalas,                                            *
        !   NIST                                                  May 2011     *
        !                                                                      *
        !***********************************************************************
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: einstein_a_au_to_si, hartree_to_wavenumber_cm, vacuum_speed_of_light_au
        use transition_output_cutoff, only: transition_probability_cutoff
        use transition_matrix_workspace, only: transition_multipole_rank, transition_multipole_type, transition_option_flags
        use level_symmetry_metadata, only: final_level_two_j_plus_one, initial_level_two_j_plus_one
        use transition_lsj_labels, only: rlev_ener_1, rlev_ener_2, string_csf1, string_csf2
        implicit none (type, external)

        ! Arguments
        integer  :: i, j, inum_ii, inum_ff
        real(real64), intent(in) :: coulomb_or_magnetic_amplitude, babushkin_amplitude
        real(real64), intent(in) :: transition_angular_frequency, factor

        ! Local declarations
        !      INTEGER :: IPAR, JPAR, ISIGNA, ISIGNB
        character (len=1) :: im
        real(real64) :: dblfac(10)
        real(real64) :: d,dd,ab,angs,angsa,ac,am,ams,sa,sb,oscc,oscb,  &
            sigma,faau,stfac,dll1,omc,fosc,eng,ggfactor,   &
            acsq,absq,oscm

        ! Initial values
        !
        !  *****  DBLFAC(I) = (2I+1)!!
        !
        data dblfac/ 3.0000000000_real64, 1.5000000000e1_real64, 1.0500000000e2_real64, &
            9.4500000000e2_real64, 1.0395000000e4_real64, 1.3513500000e5_real64, 2.0270250000e6_real64, &
            3.4459425000e7_real64, 6.5472907500e8_real64, 1.3749310575e10_real64/
        !-----------------------------------------------
        !
        d    = dabs(rlev_ener_2(inum_ff) - rlev_ener_1(inum_ii))
        if(dabs(d - dabs(transition_angular_frequency)) > 0.00000001_real64) then
            print*, "                                                      "
            print*, " INCORRECT INPUT, PROGRAM STOP!!!                     "
            print*, "                                                      "
            print*, " lsj.lbl files from the jj2lsj runs are inconsistent  "
            print*, " with the files used in the bioscl calculation.       "
            print*, " The bioscl calculation was for rci wave functions    "
            print*, " but the jj2lsj run was for rscf wave functions or    "
            print*, " vice versa.                                          "
            print*, "                                                      "
            print*, " If the bioscl calculation was for RCI wave functions "
            print*, " run jj2lsj to produce labels for RCI wave functions. "
            print*, " Then rerun rtransition.                              "
            print*, "                                                      "
            print*, " If the bioscl calculation was for rscf wave functions"
            print*, " run jj2lsj to produce labels for rscf wave functions."
            print*, " Then rerun rtransition.                              "
            stop
        end if
        dd    = d*hartree_to_wavenumber_cm
        angs  =  1.0e8_real64 / dd
        angsa = angs
        if(angs > 2000.0_real64) then
            sigma = (1.0e8_real64/angs)**2
            angsa = angs/(1.0_real64 + 8342.13e-8_real64 + 206030.0_real32/(130.0e8_real64 - sigma) &
                +15997.0_real32/(38.9e8_real64 - sigma))
        end if
        !
        !   Evaluate statistical factors and constants
        !
        stfac = initial_level_two_j_plus_one(i)
        dll1 = dble (transition_multipole_rank+transition_multipole_rank+1)
        stfac = stfac/dll1
        !
        omc = transition_angular_frequency/vacuum_speed_of_light_au
        faau = 2.0_real64*omc*stfac/dble (final_level_two_j_plus_one(j))
        fosc = vacuum_speed_of_light_au*stfac/omc
        eng = transition_angular_frequency*factor
        if (transition_option_flags(1)) eng = factor/transition_angular_frequency
        !
        ggfactor = vacuum_speed_of_light_au**(2*transition_multipole_rank-2)*dble(transition_multipole_rank)* &
            dblfac(transition_multipole_rank)**2/    &
            ((2.0_real64 * dble(transition_multipole_rank)+1.0_real64)*            &
            (dble(transition_multipole_rank)+1.0_real64)*abs(transition_angular_frequency)**(2*transition_multipole_rank-1))
        !
        !   Calculate Einstein A and B coefficients and oscillator strengths
        !
        if (transition_multipole_type == 0) then
            !
            !   Electric multipoles
            !
            !   In atomic units
            !
            acsq = coulomb_or_magnetic_amplitude**2
            absq = babushkin_amplitude**2
            ac = acsq*faau
            ab = absq*faau
            oscc = acsq*fosc
            oscb = absq*fosc
            sa=oscc*ggfactor
            sb=oscb*ggfactor
            !
            !   Convert to SI units if option 5 not set
            !
            if (.not. transition_option_flags(7)) then
                ac = ac*einstein_a_au_to_si
                ab = ab*einstein_a_au_to_si
            end if
            !
            !   Print information if both AC and AB are greater than CUTOFF
            !
            if((abs(ac)>=transition_probability_cutoff) .and. (abs(ab)>=transition_probability_cutoff)) then
                im = "E"
                if (-1.0e9_real64 < eng .and. eng < -1.0e-2_real64) then
                    write (lsj_transition_output_unit,"(A1)")  ""
                    write (lsj_transition_output_unit,"(A1)")  ""
                    if(dabs(rlev_ener_1(inum_ii)) < 9999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.8,2X,A)")                      &
                            initial_level_two_j_plus_one(i)-1,rlev_ener_1(inum_ii),              &
                            trim(string_csf1(inum_ii))
                    else if(dabs(rlev_ener_1(inum_ii)) < 99999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.7,2X,A)")                     &
                            initial_level_two_j_plus_one(i)-1,rlev_ener_1(inum_ii),             &
                            trim(string_csf1(inum_ii))
                    else
                        write(lsj_transition_output_unit,"(I4,F14.6,2X,A)")                     &
                            initial_level_two_j_plus_one(i)-1,rlev_ener_1(inum_ii),             &
                            trim(string_csf1(inum_ii))
                    end if
                    if(dabs(rlev_ener_2(inum_ff)) < 9999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.8,2X,A)")                      &
                            final_level_two_j_plus_one(j)-1,rlev_ener_2(inum_ff),              &
                            trim(string_csf2(inum_ff))
                    else if(dabs(rlev_ener_2(inum_ff)) < 99999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.7,2X,A)")                      &
                            final_level_two_j_plus_one(j)-1,rlev_ener_2(inum_ff),              &
                            trim(string_csf2(inum_ff))
                    else
                        write(lsj_transition_output_unit,"(I4,F14.6,2X,A)")                      &
                            final_level_two_j_plus_one(j)-1,rlev_ener_2(inum_ff),              &
                            trim(string_csf2(inum_ff))
                    end if
                    write (lsj_transition_output_unit,'(f11.2," CM-1",2x,f11.2," ANGS(VAC)",2x,f11.2," ANGS(AIR)"/1x,a1,i1,2x,&
                        &"S = ",1pd12.5,3x,&
                        &"GF = ",d12.5,3x,"AKI = ",d12.5,3x,"dT = ",0pf8.5)') dd,angs,angsa,                        &
                        im,transition_multipole_rank,-sb,-oscb,-ab,                            &
                        dabs(dabs(sb)-dabs(sa))/max(dabs(sb),dabs(sa))
                    write (lsj_transition_output_unit,'(9x,1pd12.5,8x,d12.5,9x,d12.5)')  -sa,-oscc,-ac
                else if ( 1.0e9_real64 > eng .and. eng > 1.0e-2_real64) then
                    write (lsj_transition_output_unit,"(A1)")  ""
                    write (lsj_transition_output_unit,"(A1)")  ""
                    if(dabs(rlev_ener_2(inum_ff)) < 9999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.8,2X,A)")                      &
                            final_level_two_j_plus_one(j)-1,rlev_ener_2(inum_ff),              &
                            trim(string_csf2(inum_ff))
                    else if(dabs(rlev_ener_2(inum_ff)) < 99999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.7,2X,A)")                      &
                            final_level_two_j_plus_one(j)-1,rlev_ener_2(inum_ff),              &
                            trim(string_csf2(inum_ff))
                    else
                        write(lsj_transition_output_unit,"(I4,F14.6,2X,A)")                      &
                            final_level_two_j_plus_one(j)-1,rlev_ener_2(inum_ff),              &
                            trim(string_csf2(inum_ff))
                    end if
                    if(dabs(rlev_ener_1(inum_ii)) < 9999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.8,2X,A)")                      &
                            initial_level_two_j_plus_one(i)-1,rlev_ener_1(inum_ii),              &
                            trim(string_csf1(inum_ii))
                    else if(dabs(rlev_ener_1(inum_ii)) < 99999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.7,2X,A)")                     &
                            initial_level_two_j_plus_one(i)-1,rlev_ener_1(inum_ii),             &
                            trim(string_csf1(inum_ii))
                    else
                        write(lsj_transition_output_unit,"(I4,F14.6,2X,A)")                     &
                            initial_level_two_j_plus_one(i)-1,rlev_ener_1(inum_ii),             &
                            trim(string_csf1(inum_ii))
                    end if
                    write (lsj_transition_output_unit,'(f11.2," CM-1",2x,f11.2," ANGS(VAC)",2x,f11.2," ANGS(AIR)"/1x,a1,i1,2x,&
                        &"S = ",1pd12.5,3x,&
                        &"GF = ",d12.5,3x,"AKI = ",d12.5,3x,"dT = ",0pf8.5)') dd,angs,angsa,                        &
                        im,transition_multipole_rank,sb,oscb,ab* &
                            final_level_two_j_plus_one(j)/initial_level_two_j_plus_one(i),       &
                        dabs(dabs(sb)-dabs(sa))/max(dabs(sb),dabs(sa))
                    write (lsj_transition_output_unit,'(9x,1pd12.5,8x,d12.5,9x,d12.5)')  sa,oscc,ac* &
                        final_level_two_j_plus_one(j)/initial_level_two_j_plus_one(i)
                end if
            end if
        else
            !
            !   Magnetic multipoles
            !
            !   In atomic units
            !
            ams = coulomb_or_magnetic_amplitude**2
            am = ams*faau
            oscm = ams*fosc
            sa=oscm*ggfactor*vacuum_speed_of_light_au*vacuum_speed_of_light_au*4
            !
            !   Convert to SI units if option 5 not set
            !
            if (.not. transition_option_flags(7)) am = am*einstein_a_au_to_si
            !
            !   Print information if AM is greater than CUTOFF
            !
            if (abs(am) >= transition_probability_cutoff) then
                im = "M"
                if (-1.0e9_real64 < eng .and. eng < -1.0e-2_real64) then
                    write (lsj_transition_output_unit,"(A1)")  ""
                    write (lsj_transition_output_unit,"(A1)")  ""
                    if(dabs(rlev_ener_1(inum_ii)) < 9999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.8,2X,A)")                      &
                            initial_level_two_j_plus_one(i)-1,rlev_ener_1(inum_ii),              &
                            trim(string_csf1(inum_ii))
                    else if(dabs(rlev_ener_1(inum_ii)) < 99999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.7,2X,A)")                     &
                            initial_level_two_j_plus_one(i)-1,rlev_ener_1(inum_ii),             &
                            trim(string_csf1(inum_ii))
                    else
                        write(lsj_transition_output_unit,"(I4,F14.6,2X,A)")                     &
                            initial_level_two_j_plus_one(i)-1,rlev_ener_1(inum_ii),             &
                            trim(string_csf1(inum_ii))
                    end if
                    if(dabs(rlev_ener_2(inum_ff)) < 9999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.8,2X,A)")                      &
                            final_level_two_j_plus_one(j)-1,rlev_ener_2(inum_ff),              &
                            trim(string_csf2(inum_ff))
                    else if(dabs(rlev_ener_2(inum_ff)) < 99999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.7,2X,A)")                      &
                            final_level_two_j_plus_one(j)-1,rlev_ener_2(inum_ff),              &
                            trim(string_csf2(inum_ff))
                    else
                        write(lsj_transition_output_unit,"(I4,F14.6,2X,A)")                      &
                            final_level_two_j_plus_one(j)-1,rlev_ener_2(inum_ff),              &
                            trim(string_csf2(inum_ff))
                    end if
                    write (lsj_transition_output_unit,'(f11.2," CM-1",2x,f11.2," ANGS(VAC)",2x,f11.2," ANGS(AIR)"/1x,a1,i1,2x,&
                        &"S = ",1pd12.5,3x,&
                        &"GF = ",d12.5,3x,"AKI = ",d12.5)') dd,angs,angsa,im,transition_multipole_rank,-sa,-oscm,-am
                else if ( 1.0e9_real64 > eng .and. eng > 1.0e-2_real64) then
                    write (lsj_transition_output_unit,"(A1)")  ""
                    write (lsj_transition_output_unit,"(A1)")  ""
                    if(dabs(rlev_ener_2(inum_ff)) < 9999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.8,2X,A)")                      &
                            final_level_two_j_plus_one(j)-1,rlev_ener_2(inum_ff),              &
                            trim(string_csf2(inum_ff))
                    else if(dabs(rlev_ener_2(inum_ff)) < 99999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.7,2X,A)")                      &
                            final_level_two_j_plus_one(j)-1,rlev_ener_2(inum_ff),              &
                            trim(string_csf2(inum_ff))
                    else
                        write(lsj_transition_output_unit,"(I4,F14.6,2X,A)")                      &
                            final_level_two_j_plus_one(j)-1,rlev_ener_2(inum_ff),              &
                            trim(string_csf2(inum_ff))
                    end if
                    if(dabs(rlev_ener_1(inum_ii)) < 9999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.8,2X,A)")                      &
                            initial_level_two_j_plus_one(i)-1,rlev_ener_1(inum_ii),              &
                            trim(string_csf1(inum_ii))
                    else if(dabs(rlev_ener_1(inum_ii)) < 99999.5_real32) then
                        write(lsj_transition_output_unit,"(I4,F14.7,2X,A)")                     &
                            initial_level_two_j_plus_one(i)-1,rlev_ener_1(inum_ii),             &
                            trim(string_csf1(inum_ii))
                    else
                        write(lsj_transition_output_unit,"(I4,F14.6,2X,A)")                     &
                            initial_level_two_j_plus_one(i)-1,rlev_ener_1(inum_ii),             &
                            trim(string_csf1(inum_ii))
                    end if
                    write (lsj_transition_output_unit,'(f11.2," CM-1",2x,f11.2," ANGS(VAC)",2x,f11.2," ANGS(AIR)"/1x,a1,i1,2x,&
                        &"S = ",1pd12.5,3x,&
                        &"GF = ",d12.5,3x,"AKI = ",d12.5)') dd,angs,angsa,                        &
                        im,transition_multipole_rank,sa,oscm,am*final_level_two_j_plus_one(j)/initial_level_two_j_plus_one(i)
                end if
            end if
        end if
        return
    end subroutine write_lsj_transition_probability
end module write_lsj_transition_probability_interface
