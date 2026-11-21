module calculate_isotope_shifts_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_isotope_shifts
contains

    !***********************************************************************
    !                                                                      *
    subroutine calculate_isotope_shifts (name)
        ! Module dependencies
        use calculation_file_units, only: one_body_cache_unit, summary_output_unit, two_body_cache_unit
        !                                                                      *
        !   This routine controls the main sequence of routine calls for the   *
        !   calculation  of the MS parameters, the electron density at the     *
        !   origin and radial expectation values.                              *
        !                                                                      *
        !   Call(s) to: [LIB92]: ALLOC, DALLOC, GETYN                          *
        !               [SMS92]: RINTDENS, VINTI                               *
        !                                                                      *
        !   Written by Per Jonsson                                             *
        !                                                                      *
        !                                         Last revision: 10 Nov 1995   *
        !                                                                      *
        !   Modified by C. Naz\'e  Feb. 2011                                   *
        !   Modified by J. Ekman   Jan. 2014                                   *
        !                                                                      *
        !***********************************************************************
        !...Translated by Gediminas Gaigalas 11/18/19
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: max_orbitals
        use array_memory_management, only: alloc, dalloc
        use debug_output_settings, only: matrix_element_screening_cutoff
        use atomic_calculation_state, only: electron_mass_u, hartree_to_wavenumber_cm, max_radial_points, nuclear_charge, &
            speed_of_light_cm_per_s
        use radial_grid_state, only: radial_grid
        use nuclear_model_parameters_state, only: nuclear_model_parameters
        use selected_level_indices, only: active_level_indices, active_level_count
        use level_symmetry_metadata, only: level_parity_signs, level_two_j_plus_one
        use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_kappa
        use isotope_shift_options, only: contact_density_by_level, total_nms_by_level, radius_by_level, &
            & radius_squared_by_level, inverse_radius_by_level, inverse_radius_squared_by_level, leading_nms_by_level, &
            & leading_sms_by_level, total_sms_by_level
        use angular_symmetry_labels, only: format_j_label, parity_labels
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues

        ! Procedure interfaces
        use field_shift_radial_integral_interface, only: field_shift_radial_integral

        ! Module dependencies
        use one_electron_radial_integral, only: integrate_one_electron_radial
        use weighted_radial_overlap, only: integrate_weighted_radial_overlap

        ! Procedure interfaces
        use field_shift_radial_integral_vector_interface, only: field_shift_radial_integral_vector
        use normal_mass_shift_radial_integral_interface, only: normal_mass_shift_radial_integral

        ! Module dependencies
        use reduced_angular_matrix_element, only: reduced_angular_element

        ! Procedure interfaces
        use vinti_radial_integral_interface, only: vinti_radial_integral
        use specific_mass_shift_radial_integral_2_interface, only: specific_mass_shift_radial_integral_2
        use specific_mass_shift_radial_integral_3_interface, only: specific_mass_shift_radial_integral_3
        use open_isotope_shift_angular_cache_interface, only: open_isotope_shift_angular_cache

        ! Module dependencies
        use input_error_handling, only: fail_input_io
        use yes_no_input, only: read_yes_no

        ! Procedure interfaces
        use read_field_shift_factors_interface, only: read_field_shift_factors
        use calculate_field_shift_factors_interface, only: calculate_field_shift_factors
        use read_specific_mass_shift_interface, only: read_specific_mass_shift
        use calculate_specific_mass_shift_interface, only: calculate_specific_mass_shift
        use fit_electron_density_interface, only: fit_electron_density
        use calculate_seltzer_field_shift_factors_interface, only: calculate_seltzer_field_shift_factors
        use read_seltzer_field_shift_factors_interface, only: read_seltzer_field_shift_factors

        ! Module dependencies
        use run_timer_start, only: begin_run_timer
        use run_timer_finish, only: finish_run_timer
        implicit none (type, external)

        ! Arguments
        character(len=24), intent(in) :: name

        ! Local declarations
        integer :: first_name_blank
        !      REAL(DOUBLE), INTENT(IN) :: DR2
        !      INTEGER, INTENT(IN)      :: NOPAR
        real(real64), pointer :: field_shift_radial_integral_profiles(:, :, :)
        real(real64), pointer   :: level_electron_density_profiles(:, :)
        real(real64), pointer     :: electron_density_fit(:)
        real(real64), pointer     :: density_fit_parameters(:, :)
        real(real64), pointer     :: effective_electron_density(:)
        real(real64), pointer     :: density_fit_residual(:)


        real(real64) :: vinti_integrals(max_orbitals, max_orbitals), recoil_corrected_vinti_integrals(max_orbitals, max_orbitals)
        real(real64) :: contact_density_integrals(max_orbitals, max_orbitals), total_nms_integrals(max_orbitals, &
            & max_orbitals), radius_integrals(max_orbitals, max_orbitals), &
            radius_squared_integrals(max_orbitals, max_orbitals), inverse_radius_integrals(max_orbitals, max_orbitals), &
            inverse_radius_squared_integrals(max_orbitals, max_orbitals), &
            leading_nms_integrals(max_orbitals, max_orbitals)
        real(real64) :: reduced_rank1_angular_element, specific_mass_shift_hz_u, normal_mass_shift_hz_u, &
            nuclear_density_fit_radius
        logical :: yes, yes2, avail_tb, avail_ob
        integer :: i, orbital_index, start_clock_count, one_body_cache_mode, two_body_cache_mode, &
            nuclear_density_grid_point_count, cache_status
        !-----------------------------------------------
        !

        !
        !   Matrix elements smaller than CUTOFF are not accumulated
        !
        !      PARAMETER (CUTOFF = 1.0D-10)

        write(*,*) "-------------------------------"
        write(*,*) "calculate_isotope_shifts: Execution Begins ..."
        write(*,*) "-------------------------------"

        !     Determine point end point NRNUC for density fit
        nuclear_density_fit_radius = nuclear_model_parameters(1)+2.0_real64*log(3.0_real64)*nuclear_model_parameters(2)
        !PARM(1) - Fermi parameter c
        !PARM(2) - Fermi parameter a
        do i=1,max_radial_points
            if(radial_grid(i+1)>nuclear_density_fit_radius.and.radial_grid(i)<=nuclear_density_fit_radius) then
                nuclear_density_grid_point_count = i+1
                exit
            end if
        end do


        !   Allocate memory

        allocate( level_electron_density_profiles(active_level_count,nuclear_density_grid_point_count), &
            & field_shift_radial_integral_profiles(max_orbitals,max_orbitals,nuclear_density_grid_point_count) )    ! JE ADD
        allocate( electron_density_fit(nuclear_density_grid_point_count) )                                   ! JE ADD

        field_shift_radial_integral_profiles(:,:,:) = 0.0_real64                                      ! JE ADD
        level_electron_density_profiles(:,:) = 0.0_real64                                        ! JE ADD

        allocate( density_fit_parameters(6,active_level_count),effective_electron_density(active_level_count), &
            density_fit_residual(active_level_count) &
            & )                ! JE ADD

        !
        !   Allocate storage for local arrays
        !
        call alloc (leading_sms_by_level, active_level_count,"leading_sms_by_level", "calculate_isotope_shifts")
        call alloc (total_sms_by_level, active_level_count,"total_sms_by_level", "calculate_isotope_shifts")
        call alloc (contact_density_by_level,active_level_count,"contact_density_by_level","calculate_isotope_shifts")
        call alloc (total_nms_by_level,active_level_count,"total_nms_by_level","calculate_isotope_shifts")
        call alloc (radius_by_level,active_level_count,"radius_by_level","calculate_isotope_shifts")
        call alloc (radius_squared_by_level,active_level_count,"radius_squared_by_level","calculate_isotope_shifts")
        call alloc (inverse_radius_by_level,active_level_count,"inverse_radius_by_level","calculate_isotope_shifts")
        call alloc (inverse_radius_squared_by_level,active_level_count,"inverse_radius_squared_by_level","calculate_isotope_shifts")
        call alloc (leading_nms_by_level,active_level_count,"leading_nms_by_level","calculate_isotope_shifts")

        call begin_run_timer (start_clock_count, "RIS_CAL")

        !
        !   Initialise
        !
        do i = 1,active_level_count
            leading_sms_by_level(i)  = 0.0_real64
            total_sms_by_level(i)  = 0.0_real64
            contact_density_by_level(i) = 0.0_real64
            total_nms_by_level(i) = 0.0_real64
            radius_by_level(i) = 0.0_real64
            radius_squared_by_level(i) = 0.0_real64
            inverse_radius_by_level(i) = 0.0_real64
            inverse_radius_squared_by_level(i) = 0.0_real64
            leading_nms_by_level(i) = 0.0_real64
        end do

        !      VSH   = .TRUE.
        !      SMSSH = .TRUE.
        !
        !   Calculate all integrals needed for the volume shift calc.
        !
        !      IF (VSH) THEN

        write(*,*) " Compute higher order field shift electronic factors?"
        yes2 = read_yes_no ()

        do i = 1,orbital_count
            do orbital_index = 1,orbital_count
                if (orbital_kappa(i)==orbital_kappa(orbital_index)) then
                    contact_density_integrals(i,orbital_index) = field_shift_radial_integral(i,orbital_index)
                    !GG                 DINT2(I,J) = RINTI(I,J,1)
                    if(yes2) then
                        call field_shift_radial_integral_vector(i,orbital_index,field_shift_radial_integral_profiles, &
                            & nuclear_density_grid_point_count)
                    end if
                    call normal_mass_shift_radial_integral(i,orbital_index,total_nms_integrals(i,orbital_index), &
                        & leading_nms_integrals(i,orbital_index))
                    radius_integrals(i,orbital_index) = integrate_weighted_radial_overlap(i,orbital_index,1)
                    radius_squared_integrals(i,orbital_index) = integrate_weighted_radial_overlap(i,orbital_index,2)
                    inverse_radius_integrals(i,orbital_index) = integrate_weighted_radial_overlap(i,orbital_index,-1)
                    inverse_radius_squared_integrals(i,orbital_index) = integrate_weighted_radial_overlap(i,orbital_index,-2)
                else
                    contact_density_integrals(i,orbital_index) = 0.0_real64
                    field_shift_radial_integral_profiles(i,orbital_index,:) = 0.0_real64
                    total_nms_integrals(i,orbital_index) = 0.0_real64
                    radius_integrals(i,orbital_index) = 0.0_real64
                    radius_squared_integrals(i,orbital_index) = 0.0_real64
                    inverse_radius_integrals(i,orbital_index) = 0.0_real64
                    inverse_radius_squared_integrals(i,orbital_index) = 0.0_real64
                    leading_nms_integrals(i,orbital_index) = 0.0_real64
                end if
            end do
        end do
        !      ENDIF
        !
        !   Calculate and save the Vinti integrals
        !
        !      IF (SMSSH) THEN
        do i = 1,orbital_count
            do orbital_index = 1,orbital_count
                if (i/=orbital_index) then
                    reduced_rank1_angular_element = reduced_angular_element(orbital_kappa(i),1,orbital_kappa(orbital_index))
                    if (dabs(reduced_rank1_angular_element) > matrix_element_screening_cutoff) then
                        vinti_integrals  (i,orbital_index) = vinti_radial_integral(i,orbital_index)
                        recoil_corrected_vinti_integrals(i,orbital_index) = vinti_integrals(i, &
                            & orbital_index)                                   &
                            + specific_mass_shift_radial_integral_2(i,orbital_index)/ &
                                & reduced_rank1_angular_element                         &
                            + specific_mass_shift_radial_integral_3(i,orbital_index)
                    else
                        vinti_integrals (i,orbital_index) = 0.0_real64
                        recoil_corrected_vinti_integrals(i,orbital_index) = 0.0_real64
                    end if
                else
                    vinti_integrals (i,orbital_index) = 0.0_real64
                    recoil_corrected_vinti_integrals(i,orbital_index) = 0.0_real64
                end if
            end do
        end do
        !      ENDIF
        !
        !   See if the appropriate angular data is available. If so,
        !   then read the angular files and perform the calculation.
        !   If not, the user can choose to save them or not
        !
        one_body_cache_mode = 0
        two_body_cache_mode = 0
        call open_isotope_shift_angular_cache(name,avail_ob,1,cache_status,one_body_cache_unit)
        if (cache_status /= 0) call fail_input_io(cache_status,one_body_cache_unit,"ris4: one-body angular cache "//trim(name))
        call open_isotope_shift_angular_cache(name,avail_tb,2,cache_status,two_body_cache_unit)
        if (cache_status /= 0) call fail_input_io(cache_status,two_body_cache_unit,"ris4: two-body angular cache "//trim(name))
        if ((.not. avail_ob) .and. (.not. avail_tb)) then
            print *," Save ang. coefficients of one- and two-body op.?"
            yes = read_yes_no ()
            print *
            !C        IF(.NOT. YES) THEN
            !C          PRINT *,' Save ang. coefficients of one-body op. ONLY?'
            !C          YES = GETYN ()
            !C          PRINT *
            !C          IF(.NOT. YES) THEN
            !C            PRINT *,' Save ang. coefficients of two-body op. ONLY?'
            !C            YES = GETYN ()
            !C            PRINT *
            !C            IF(YES) DOIT_TB = 1
            !C          ELSE
            !C            DOIT_OB = 1
            !C          ENDIF
            !C        ELSE
            if(yes) then
                one_body_cache_mode = 1
                two_body_cache_mode = 1
            end if
        else if (.not. avail_ob) then
            print *," Save ang. coefficients of one-body op. ?"
            yes = read_yes_no ()
            print *
            if(yes) one_body_cache_mode = 1
        else if (.not. avail_tb) then
            print *," Save ang. coefficients of two-body op. ?"
            yes = read_yes_no ()
            print *
            if(yes) two_body_cache_mode = 1
        end if
        if (avail_ob) then
            !    angdata leaves each available cache connected and rewound.
            if(yes2) then
                call read_seltzer_field_shift_factors(contact_density_integrals,total_nms_integrals,radius_integrals, & ! JE ADD
                    radius_squared_integrals,inverse_radius_integrals,inverse_radius_squared_integrals, & ! JE ADD
                    leading_nms_integrals,field_shift_radial_integral_profiles,level_electron_density_profiles, &
                        & nuclear_density_grid_point_count) ! JE ADD
            else
                call read_field_shift_factors(contact_density_integrals,total_nms_integrals,radius_integrals, &
                    & radius_squared_integrals,inverse_radius_integrals,inverse_radius_squared_integrals,leading_nms_integrals)
            end if
            close(one_body_cache_unit)

        !        CALL DENSREAD(DINT1,DINT2,DINT3,             &   ! JE ADD
        !                      DINT4,DINT5,DINT6,             &   ! JE ADD
        !                      DINT7,DINT1VEC,DENS1VEC)           ! JE ADD
        else
            !Check if the user wants to save one-body ang. coeff.
            if(one_body_cache_mode ==1) then
                first_name_blank = index(name," ")
                open(newunit=one_body_cache_unit,file = name(1:first_name_blank-1)//".IOB",status="UNKNOWN"     &
                    ,form="UNFORMATTED")
            end if
            if(yes2) then
                call calculate_seltzer_field_shift_factors(one_body_cache_mode,contact_density_integrals, &
                    & total_nms_integrals,radius_integrals, & ! JE ADD
                    radius_squared_integrals,inverse_radius_integrals,inverse_radius_squared_integrals,leading_nms_integrals, & ! JE ADD
                    field_shift_radial_integral_profiles,level_electron_density_profiles,nuclear_density_grid_point_count) ! JE ADD
            else
                call calculate_field_shift_factors(one_body_cache_mode,contact_density_integrals,total_nms_integrals, &
                    & radius_integrals,radius_squared_integrals, &
                    inverse_radius_integrals,inverse_radius_squared_integrals,leading_nms_integrals)
            end if

            !        CALL NATORBNEW(NAME)                            ! JE ADD
            !        CALL DENSNEW(DOIT_OB,DINT1,DINT2,DINT3,     &   ! JE ADD
            !             DINT4,DINT5,DINT6,DINT7,               &   ! JE ADD
            !             DINT1VEC,DENS1VEC)                         ! JE ADD
        end if
        !
        !
        !
        if (avail_tb) then
            call read_specific_mass_shift(vinti_integrals,recoil_corrected_vinti_integrals)
            close(two_body_cache_unit)
        else
            ! Check if the user wants to save two-body ang. coeff.
            if(two_body_cache_mode==1) then
                first_name_blank = index(name," ")
                open(newunit=two_body_cache_unit,file = name(1:first_name_blank-1)//".ITB",status="UNKNOWN"     &
                    ,form="UNFORMATTED")
            end if
            call calculate_specific_mass_shift(two_body_cache_mode,vinti_integrals,recoil_corrected_vinti_integrals)
        end if
        !
        !   Printouts
        !

        write (summary_output_unit,"(a24,i3)") "Number of eigenvalues: ",active_level_count
        write(summary_output_unit,*)
        write(summary_output_unit,*)
        write (summary_output_unit,'(" Level  J Parity  Energy")')
        do i = 1, active_level_count
            write (summary_output_unit,'(1x,i3,5x,2a,3x,d20.10,"  (a.u.)")') active_level_indices(i), &
                format_j_label(level_two_j_plus_one(i)), &
                & parity_labels((level_parity_signs(i)+3)/2), &
                hamiltonian_energy_shift+shifted_eigenvalues(i)
        end do
        write (summary_output_unit,'(//" Level  J Parity  Normal mass shift parameter")')
        do i = 1,active_level_count
            write(summary_output_unit,'(/,29x,"<K^1>",13x,"<K^2+K^3>",9x,"<K^1+K^2+K^3>")')
            write (summary_output_unit,'(1x,i3,5x,2a,3x,3d20.10,"  (a.u.)")') active_level_indices(i), &
                format_j_label(level_two_j_plus_one(i)), &
                & parity_labels((level_parity_signs(i)+3)/2), &
                leading_nms_by_level(i),(total_nms_by_level(i)-leading_nms_by_level(i)),total_nms_by_level(i)
            normal_mass_shift_hz_u=leading_nms_by_level(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s
            !         if (dabs(HzNMSu) .LE. (10**8)) THEN
            !             WRITE (24,315) HzNMSu*(1.0D-06),
            !     :      (DENS2(I)-DENS7(I))*AUMAMU*AUCM*CCMS*(1.0D-06),
            !     :       DENS2(I)*AUMAMU*AUCM*CCMS*(1.0D-06)
            !         else
            write (summary_output_unit,'(20x,3d20.10,2x,"(GHz u)")') normal_mass_shift_hz_u*(1.0e-9_real64), &
                (total_nms_by_level(i)-leading_nms_by_level(i))*electron_mass_u*hartree_to_wavenumber_cm* &
                    & speed_of_light_cm_per_s*(1.0e-9_real64), &
                total_nms_by_level(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64)
            !         endif
        end do

        write (summary_output_unit,'(//" Level  J Parity  Specific mass shift parameter")')
        do i = 1,active_level_count
            write(summary_output_unit,'(/,29x,"<K^1>",13x,"<K^2+K^3>",9x,"<K^1+K^2+K^3>")')
            write (summary_output_unit,'(1x,i3,5x,2a,3x,3d20.10,"  (a.u.)")') active_level_indices(i), &
                format_j_label(level_two_j_plus_one(i)), &
                & parity_labels((level_parity_signs(i)+3)/2), &
                leading_sms_by_level(i),total_sms_by_level(i)-leading_sms_by_level(i),total_sms_by_level(i)
            specific_mass_shift_hz_u=leading_sms_by_level(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s
            if (dabs(specific_mass_shift_hz_u) <= (10**8)) then
                write (summary_output_unit,'(20x,3d20.10,2x,"(MHz u)")') specific_mass_shift_hz_u*(1.0e-6_real64), &
                    (total_sms_by_level(i)-leading_sms_by_level(i))*electron_mass_u*hartree_to_wavenumber_cm* &
                        & speed_of_light_cm_per_s*(1.0e-6_real64), &
                    total_sms_by_level(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-6_real64)
            else
                write (summary_output_unit,'(20x,3d20.10,2x,"(GHz u)")') specific_mass_shift_hz_u*(1.0e-9_real64), &
                    (total_sms_by_level(i)-leading_sms_by_level(i))*electron_mass_u*hartree_to_wavenumber_cm* &
                        & speed_of_light_cm_per_s*(1.0e-9_real64), &
                    total_sms_by_level(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64)
            end if
        end do
        if(yes2) then
            write (summary_output_unit,'(//" Level  J Parity  Electron density in atomic units"//24x,"Dens. (a.u.)")')
            do i = 1,active_level_count
                call fit_electron_density(radial_grid,level_electron_density_profiles(i,:),nuclear_charge, &
                    & nuclear_model_parameters,nuclear_density_grid_point_count,density_fit_parameters(:,i), &
                    & effective_electron_density(i), &
                    & density_fit_residual(i))  ! JE ADD
                write (summary_output_unit,'(1x,i3,5x,2a,3x,d20.10)') active_level_indices(i), &
                    format_j_label(level_two_j_plus_one(i)), &
                    parity_labels((level_parity_signs(i)+3)/2),effective_electron_density(i)
            end do
            write (summary_output_unit, &
                & '(//" Level  J Parity  Field shift electronic factors and average point discrepancy in fit"//&
                &24x,"F0 (GHz/fm^2)",6x," F2 (GHz/fm^4)",6x," F4 (GHz/fm^6)",6x," F6 (GHz/fm^8)",6x,&
                &" Disc. (per mille)")')
            do i = 1,active_level_count
                !CALL EDENSITYFIT(R,DENS1VEC(I,:),Z,PARM,NRNUC,FMAT,RHO,RES)  ! JE ADD
                write (summary_output_unit,'(1x,i3,5x,2a,3x,4d20.10,f10.4)') active_level_indices(i), &
                    format_j_label(level_two_j_plus_one(i)), &
                    parity_labels((level_parity_signs(i)+3)/2),density_fit_parameters(1,i),density_fit_parameters(2,i), &
                    density_fit_parameters(3,i),density_fit_parameters(4,i),density_fit_residual(i)
            end do
            write (summary_output_unit,'(//&
                &" Level  J Parity  Field shift electronic factors (corrected for varying density inside nucleus)"//&
                &24x,"F0VED0 (GHz/fm^2)   F0VED1 (GHz/fm^4)")')
            do i = 1,active_level_count
                !CALL EDENSITYFIT(R,DENS1VEC(I,:),Z,PARM,NRNUC,FMAT,RHO,RES)  ! JE ADD
                write (summary_output_unit,'(1x,i3,5x,2a,3x,2d20.10)') active_level_indices(i), &
                    format_j_label(level_two_j_plus_one(i)), &
                    parity_labels((level_parity_signs(i)+3)/2),density_fit_parameters(5,i),density_fit_parameters(6,i)
            end do
            write(summary_output_unit,*)
        else
            write (summary_output_unit,'(//" Electron density in atomic units"//" Level  J Parity",8x,"DENS (a.u.)"/)')
            do i = 1,active_level_count
                write (summary_output_unit,'(1x,i3,5x,2a,3x,d20.10)') active_level_indices(i), &
                    format_j_label(level_two_j_plus_one(i)), &
                    parity_labels((level_parity_signs(i)+3)/2),contact_density_by_level(i)
            end do
        end if
        !     WRITE (24,309)
        !      DO 20 I = 1,NVEC
        !     WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
        !     :                  DENS3(I)
        !     20 CONTINUE
        !     WRITE (24,310)
        !     DO 22 I = 1,NVEC
        !     WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
        !     :                  DENS4(I)
        !     22 CONTINUE
        !     WRITE (24,311)
        !     DO 24 I = 1,NVEC
        !         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
        !     :                  DENS5(I)
        !     24 CONTINUE
        !     WRITE (24,312)
        !     DO 26 I = 1,NVEC
        !     WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
        !     :                  DENS6(I)
        !     26 CONTINUE
        !
        !   Dealloc
        !
        call dalloc (leading_sms_by_level,"leading_sms_by_level","calculate_isotope_shifts")
        call dalloc (total_sms_by_level,"total_sms_by_level","calculate_isotope_shifts")
        call dalloc (contact_density_by_level,"contact_density_by_level","calculate_isotope_shifts")
        call dalloc (total_nms_by_level,"total_nms_by_level","calculate_isotope_shifts")
        call dalloc (radius_by_level,"radius_by_level","calculate_isotope_shifts")
        call dalloc (radius_squared_by_level,"radius_squared_by_level","calculate_isotope_shifts")
        call dalloc (inverse_radius_by_level,"inverse_radius_by_level","calculate_isotope_shifts")
        call dalloc (inverse_radius_squared_by_level,"inverse_radius_squared_by_level","calculate_isotope_shifts")
        call dalloc (leading_nms_by_level,"leading_nms_by_level","calculate_isotope_shifts")

        write(*,*) "-------------------------------"
        write(*,*) "calculate_isotope_shifts: Execution Finished ..."
        write(*,*) "-------------------------------"

        call finish_run_timer (start_clock_count, "RIS_CAL")
        return
        !

        !     Electron density values

        !     Electronic factor values

        !     Field shift values

        !     Electron density header


        !     Electronic factors header

        !     Electronic field shift header



        !     309 FORMAT (//' Radial expectationvalue'
        !     :        //' Level  J Parity',8X,'<r> (a.u.)'/)
        !     310 FORMAT (//' Radial expectationvalue'
        !     :        //' Level  J Parity',8X,'<r2> (a.u.)'/)
        !     311 FORMAT (//' Radial expectationvalue'
        !     :        //' Level  J Parity',8X,'<r-1> (a.u.)'/)
        !     312 FORMAT (//' Radial expectationvalue'
        !     :        //' Level  J Parity',8X,'<r-2> (a.u.)'/)
        !
        return
    end subroutine calculate_isotope_shifts
end module calculate_isotope_shifts_interface
