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
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
    use atomic_calculation_state, only: electron_mass_u, hartree_to_wavenumber_cm, max_radial_points, nuclear_charge, &
        speed_of_light_cm_per_s
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
    use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
    use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
    use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_kappa
    use isotope_shift_integral_tables, only: ndtea, ntei, pindte, pvalte, firstt
    use angular_coefficient_buffer_state, only: angular_buffer_capacity, angular_coefficient_count, &
        angular_coefficient_labels, angular_coefficients
    use isotope_shift_options, only: dens1, dens2, dens3, dens4, dens5, dens6, dens7, smsc1, smsc2
    use angular_symmetry_labels, only: right_justified_j_labels => jlbr, parity_labels => jlbp
    use isotope_shift_options, only: dens1, dens2, dens3, dens4, dens5, dens6, dens7, smsc1, smsc2
    use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
        csf_mixing_coefficients, evecff, evecii

    ! Procedure interfaces
    use field_shift_radial_integral_interface, only: field_shift_radial_integral

    ! Module dependencies
    use one_electron_radial_integral, only: rinti
    use weighted_radial_overlap, only: rint

    ! Procedure interfaces
    use field_shift_radial_integral_vector_interface, only: field_shift_radial_integral_vector
    use normal_mass_shift_radial_integral_interface, only: normal_mass_shift_radial_integral

    ! Module dependencies
    use reduced_angular_matrix_element, only: cre

    ! Procedure interfaces
    use vinti_radial_integral_interface, only: vinti_radial_integral
    use specific_mass_shift_radial_integral_2_interface, only: specific_mass_shift_radial_integral_2
    use specific_mass_shift_radial_integral_3_interface, only: specific_mass_shift_radial_integral_3
    use open_isotope_shift_angular_cache_interface, only: open_isotope_shift_angular_cache

    ! Module dependencies
    use input_error_handling, only: fail_input_io
    use yes_no_input, only: getyn

    ! Procedure interfaces
    use read_field_shift_factors_interface, only: read_field_shift_factors
    use calculate_field_shift_factors_interface, only: calculate_field_shift_factors
    use read_specific_mass_shift_interface, only: read_specific_mass_shift
    use calculate_specific_mass_shift_interface, only: calculate_specific_mass_shift
    use fit_electron_density_interface, only: fit_electron_density
    use calculate_seltzer_field_shift_factors_interface, only: calculate_seltzer_field_shift_factors
    use read_seltzer_field_shift_factors_interface, only: read_seltzer_field_shift_factors

    ! Module dependencies
    use run_timer_start, only: starttime
    use run_timer_finish, only: stoptime
    implicit none (type, external)

    ! Arguments
    character(len=24), intent(in) :: name

    ! Local declarations
    !      REAL(DOUBLE), INTENT(IN) :: DR2
    !      INTEGER, INTENT(IN)      :: NOPAR
    real(real64), pointer :: dint1vec(:, :, :)
    real(real64), pointer   :: dens1vec(:, :)
    real(real64), pointer     :: densfit(:)
    real(real64), pointer     :: fmat(:, :)
    real(real64), pointer     :: rho(:)
    real(real64), pointer     :: res(:)


    real(real64)      :: shell_coefficients(max_orbitals)
    real(real64) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
    real(real64) :: dint1(max_orbitals, max_orbitals), dint2(max_orbitals, max_orbitals), dint3(max_orbitals, max_orbitals), &
        dint4(max_orbitals, max_orbitals), dint5(max_orbitals, max_orbitals), &
        dint6(max_orbitals, max_orbitals), &
        dint7(max_orbitals, max_orbitals)
    real(real64) :: au2fm, rcre, hzsmsu, hznmsu, fo90
    character(len=11) :: number_text
    character(len=2)  :: ck
    logical :: vsh, nucde, smssh, yes, yes2, avail_tb, avail_ob
    integer :: i, orbital_index, l, start_clock_count, doit_ob, doit_tb, nrnuc, cache_status
    !-----------------------------------------------
    !

    !
    !   Matrix elements smaller than CUTOFF are not accumulated
    !
    !      PARAMETER (CUTOFF = 1.0D-10)

    write(*,*) "-------------------------------"
    write(*,*) "RIS_CAL: Execution Begins ..."
    write(*,*) "-------------------------------"

    !     Determine point end point NRNUC for density fit
    fo90 = nuclear_model_parameters(1)+2.0_real64*log(3.0_real64)*nuclear_model_parameters(2)
    !PARM(1) - Fermi parameter c
    !PARM(2) - Fermi parameter a
    do i=1,max_radial_points
        if(radial_grid(i+1)>fo90.and.radial_grid(i)<=fo90) then
            nrnuc = i+1
            exit
        end if
    end do


    !   Allocate memory

    allocate( dens1vec(nvec,nrnuc),dint1vec(max_orbitals,max_orbitals,nrnuc) )    ! JE ADD
    allocate( densfit(nrnuc) )                                   ! JE ADD

    dint1vec(:,:,:) = 0.0_real64                                      ! JE ADD
    dens1vec(:,:) = 0.0_real64                                        ! JE ADD

    allocate( fmat(6,nvec),rho(nvec),res(nvec) )                ! JE ADD

    !
    !   Allocate storage for local arrays
    !
    call alloc (smsc1, nvec,"SMSC1", "RIS_CAL")
    call alloc (smsc2, nvec,"SMSC2", "RIS_CAL")
    call alloc (dens1,nvec,"DENS1","RIS_CAL")
    call alloc (dens2,nvec,"DENS2","RIS_CAL")
    call alloc (dens3,nvec,"DENS3","RIS_CAL")
    call alloc (dens4,nvec,"DENS4","RIS_CAL")
    call alloc (dens5,nvec,"DENS5","RIS_CAL")
    call alloc (dens6,nvec,"DENS6","RIS_CAL")
    call alloc (dens7,nvec,"DENS7","RIS_CAL")

    call starttime (start_clock_count, "RIS_CAL")

    !
    !   Initialise
    !
    do i = 1,nvec
        smsc1(i)  = 0.0_real64
        smsc2(i)  = 0.0_real64
        dens1(i) = 0.0_real64
        dens2(i) = 0.0_real64
        dens3(i) = 0.0_real64
        dens4(i) = 0.0_real64
        dens5(i) = 0.0_real64
        dens6(i) = 0.0_real64
        dens7(i) = 0.0_real64
    end do

    !      VSH   = .TRUE.
    !      SMSSH = .TRUE.
    !
    !   Calculate all integrals needed for the volume shift calc.
    !
    !      IF (VSH) THEN

    write(*,*) " Compute higher order field shift electronic factors?"
    yes2 = getyn ()

    do i = 1,orbital_count
        do orbital_index = 1,orbital_count
            if (orbital_kappa(i)==orbital_kappa(orbital_index)) then
                dint1(i,orbital_index) = field_shift_radial_integral(i,orbital_index)
                !GG                 DINT2(I,J) = RINTI(I,J,1)
                if(yes2) then
                    call field_shift_radial_integral_vector(i,orbital_index,dint1vec,nrnuc)
                end if
                call normal_mass_shift_radial_integral(i,orbital_index,dint2(i,orbital_index),dint7(i,orbital_index))
                dint3(i,orbital_index) = rint(i,orbital_index,1)
                dint4(i,orbital_index) = rint(i,orbital_index,2)
                dint5(i,orbital_index) = rint(i,orbital_index,-1)
                dint6(i,orbital_index) = rint(i,orbital_index,-2)
            else
                dint1(i,orbital_index) = 0.0_real64
                dint1vec(i,orbital_index,:) = 0.0_real64
                dint2(i,orbital_index) = 0.0_real64
                dint3(i,orbital_index) = 0.0_real64
                dint4(i,orbital_index) = 0.0_real64
                dint5(i,orbital_index) = 0.0_real64
                dint6(i,orbital_index) = 0.0_real64
                dint7(i,orbital_index) = 0.0_real64
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
                rcre = cre(orbital_kappa(i),1,orbital_kappa(orbital_index))
                if (dabs(rcre) > cutoff) then
                    vint  (i,orbital_index) = vinti_radial_integral(i,orbital_index)
                    vint2(i,orbital_index) = vint(i,orbital_index)                                   &
                        + specific_mass_shift_radial_integral_2(i,orbital_index)/rcre                         &
                        + specific_mass_shift_radial_integral_3(i,orbital_index)
                else
                    vint (i,orbital_index) = 0.0_real64
                    vint2(i,orbital_index) = 0.0_real64
                end if
            else
                vint (i,orbital_index) = 0.0_real64
                vint2(i,orbital_index) = 0.0_real64
            end if
        end do
    end do
    !      ENDIF
    !
    !   See if the appropriate angular data is available. If so,
    !   then read the angular files and perform the calculation.
    !   If not, the user can choose to save them or not
    !
    doit_ob = 0
    doit_tb = 0
    call open_isotope_shift_angular_cache(name,avail_ob,1,cache_status,one_body_cache_unit)
    if (cache_status /= 0) call fail_input_io(cache_status,one_body_cache_unit,"ris4: one-body angular cache "//trim(name))
    call open_isotope_shift_angular_cache(name,avail_tb,2,cache_status,two_body_cache_unit)
    if (cache_status /= 0) call fail_input_io(cache_status,two_body_cache_unit,"ris4: two-body angular cache "//trim(name))
    if ((.not. avail_ob) .and. (.not. avail_tb)) then
        print *," Save ang. coefficients of one- and two-body op.?"
        yes = getyn ()
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
            doit_ob = 1
            doit_tb = 1
        end if
    else if (.not. avail_ob) then
        print *," Save ang. coefficients of one-body op. ?"
        yes = getyn ()
        print *
        if(yes) doit_ob = 1
    else if (.not. avail_tb) then
        print *," Save ang. coefficients of two-body op. ?"
        yes = getyn ()
        print *
        if(yes) doit_tb = 1
    end if
    if (avail_ob) then
        !    angdata leaves each available cache connected and rewound.
        if(yes2) then
            call read_seltzer_field_shift_factors(dint1,dint2,dint3, & ! JE ADD
                dint4,dint5,dint6, & ! JE ADD
                dint7,dint1vec,dens1vec,nrnuc) ! JE ADD
        else
            call read_field_shift_factors(dint1,dint2,dint3,dint4,dint5,dint6,dint7)
        end if
        close(one_body_cache_unit)

    !        CALL DENSREAD(DINT1,DINT2,DINT3,             &   ! JE ADD
    !                      DINT4,DINT5,DINT6,             &   ! JE ADD
    !                      DINT7,DINT1VEC,DENS1VEC)           ! JE ADD
    else
        !Check if the user wants to save one-body ang. coeff.
        if(doit_ob ==1) then
            orbital_index = index(name," ")
            open(newunit=one_body_cache_unit,file = name(1:orbital_index-1)//".IOB",status="UNKNOWN"     &
                ,form="UNFORMATTED")
        end if
        if(yes2) then
            call calculate_seltzer_field_shift_factors(doit_ob,dint1,dint2,dint3, & ! JE ADD
                dint4,dint5,dint6,dint7, & ! JE ADD
                dint1vec,dens1vec,nrnuc) ! JE ADD
        else
            call calculate_field_shift_factors(doit_ob,dint1,dint2,dint3,dint4, &
                dint5,dint6,dint7)
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
        call read_specific_mass_shift(vint,vint2)
        close(two_body_cache_unit)
    else
        ! Check if the user wants to save two-body ang. coeff.
        if(doit_tb==1) then
            orbital_index = index(name," ")
            open(newunit=two_body_cache_unit,file = name(1:orbital_index-1)//".ITB",status="UNKNOWN"     &
                ,form="UNFORMATTED")
        end if
        call calculate_specific_mass_shift(doit_tb,vint,vint2)
    end if
    !
    !   Printouts
    !

    write (summary_output_unit,"(a24,i3)") "Number of eigenvalues: ",nvec
    write(summary_output_unit,*)
    write(summary_output_unit,*)
    write (summary_output_unit,'(" Level  J Parity  Energy")')
    do i = 1, nvec
        write (summary_output_unit,'(1x,i3,5x,2a4,3x,d20.10,"  (a.u.)")') ivec(i),right_justified_j_labels(iatjpo(i)), &
            & parity_labels((iaspar(i)+3)/2), &
            hamiltonian_energy_shift+shifted_eigenvalues(i)
    end do
    write (summary_output_unit,'(//" Level  J Parity  Normal mass shift parameter")')
    do i = 1,nvec
        write(summary_output_unit,'(/,29x,"<K^1>",13x,"<K^2+K^3>",9x,"<K^1+K^2+K^3>")')
        write (summary_output_unit,'(1x,i3,5x,2a4,3x,3d20.10,"  (a.u.)")') ivec(i),right_justified_j_labels(iatjpo(i)), &
            & parity_labels((iaspar(i)+3)/2), &
            dens7(i),(dens2(i)-dens7(i)),dens2(i)
        hznmsu=dens7(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s
        !         if (dabs(HzNMSu) .LE. (10**8)) THEN
        !             WRITE (24,315) HzNMSu*(1.0D-06),
        !     :      (DENS2(I)-DENS7(I))*AUMAMU*AUCM*CCMS*(1.0D-06),
        !     :       DENS2(I)*AUMAMU*AUCM*CCMS*(1.0D-06)
        !         else
        write (summary_output_unit,'(20x,3d20.10,2x,"(GHz u)")') hznmsu*(1.0e-9_real64), &
            (dens2(i)-dens7(i))*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64), &
            dens2(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64)
        !         endif
    end do

    write (summary_output_unit,'(//" Level  J Parity  Specific mass shift parameter")')
    do i = 1,nvec
        write(summary_output_unit,'(/,29x,"<K^1>",13x,"<K^2+K^3>",9x,"<K^1+K^2+K^3>")')
        write (summary_output_unit,'(1x,i3,5x,2a4,3x,3d20.10,"  (a.u.)")') ivec(i),right_justified_j_labels(iatjpo(i)), &
            & parity_labels((iaspar(i)+3)/2), &
            smsc1(i),smsc2(i)-smsc1(i),smsc2(i)
        hzsmsu=smsc1(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s
        if (dabs(hzsmsu) <= (10**8)) then
            write (summary_output_unit,'(20x,3d20.10,2x,"(MHz u)")') hzsmsu*(1.0e-6_real64), &
                (smsc2(i)-smsc1(i))*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-6_real64), &
                smsc2(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-6_real64)
        else
            write (summary_output_unit,'(20x,3d20.10,2x,"(GHz u)")') hzsmsu*(1.0e-9_real64), &
                (smsc2(i)-smsc1(i))*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64), &
                smsc2(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64)
        end if
    end do
    if(yes2) then
        write (summary_output_unit,'(//" Level  J Parity  Electron density in atomic units"//24x,"Dens. (a.u.)")')
        do i = 1,nvec
            call fit_electron_density(radial_grid,dens1vec(i,:),nuclear_charge,nuclear_model_parameters,nrnuc,fmat(:,i),rho(i), &
                & res(i))  ! JE ADD
            write (summary_output_unit,'(1x,i3,5x,2a4,3x,d20.10)') ivec(i),right_justified_j_labels(iatjpo(i)), &
                parity_labels((iaspar(i)+3)/2),rho(i)
        end do
        write (summary_output_unit,'(//" Level  J Parity  Field shift electronic factors and average point discrepancy in fit"//&
            &24x,"F0 (GHz/fm^2)",6x," F2 (GHz/fm^4)",6x," F4 (GHz/fm^6)",6x," F6 (GHz/fm^8)",6x,&
            &" Disc. (per mille)")')
        do i = 1,nvec
            !CALL EDENSITYFIT(R,DENS1VEC(I,:),Z,PARM,NRNUC,FMAT,RHO,RES)  ! JE ADD
            write (summary_output_unit,'(1x,i3,5x,2a4,3x,4d20.10,f10.4)') ivec(i),right_justified_j_labels(iatjpo(i)), &
                parity_labels((iaspar(i)+3)/2),fmat(1,i),fmat(2,i), &
                fmat(3,i),fmat(4,i),res(i)
        end do
        write (summary_output_unit,'(//&
            &" Level  J Parity  Field shift electronic factors (corrected for varying density inside nucleus)"//&
            &24x,"F0VED0 (GHz/fm^2)   F0VED1 (GHz/fm^4)")')
        do i = 1,nvec
            !CALL EDENSITYFIT(R,DENS1VEC(I,:),Z,PARM,NRNUC,FMAT,RHO,RES)  ! JE ADD
            write (summary_output_unit,'(1x,i3,5x,2a4,3x,2d20.10)') ivec(i),right_justified_j_labels(iatjpo(i)), &
                parity_labels((iaspar(i)+3)/2),fmat(5,i),fmat(6,i)
        end do
        write(summary_output_unit,*)
    else
        write (summary_output_unit,'(//" Electron density in atomic units"//" Level  J Parity",8x,"DENS (a.u.)"/)')
        do i = 1,nvec
            write (summary_output_unit,'(1x,i3,5x,2a4,3x,d20.10)') ivec(i),right_justified_j_labels(iatjpo(i)), &
                parity_labels((iaspar(i)+3)/2),dens1(i)
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
    call dalloc (smsc1,"SMSC1","RIS_CAL")
    call dalloc (smsc2,"SMSC2","RIS_CAL")
    call dalloc (dens1,"DENS1","RIS_CAL")
    call dalloc (dens2,"DENS2","RIS_CAL")
    call dalloc (dens3,"DENS3","RIS_CAL")
    call dalloc (dens4,"DENS4","RIS_CAL")
    call dalloc (dens5,"DENS5","RIS_CAL")
    call dalloc (dens6,"DENS6","RIS_CAL")
    call dalloc (dens7,"DENS7","RIS_CAL")

    write(*,*) "-------------------------------"
    write(*,*) "RIS_CAL: Execution Finished ..."
    write(*,*) "-------------------------------"

    call stoptime (start_clock_count, "RIS_CAL")
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
