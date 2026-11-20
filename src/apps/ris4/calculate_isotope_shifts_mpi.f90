!***********************************************************************
!                                                                      *
subroutine calculate_isotope_shifts_mpi (name)
    ! Module dependencies
    use calculation_file_units, only: one_body_cache_unit, summary_output_unit, two_body_cache_unit
    !                                                                      *
    !   MPI version of RIS_CAL. Each rank accumulates a strided subset of  *
    !   configuration-pair contributions and the partial sums are reduced   *
    !   before the root rank writes the final output.                       *
    !                                                                      *
    !***********************************************************************
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
    use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
        csf_mixing_coefficients, evecff, evecii
    use mpi_state, only: mpi_barrier, mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_rank

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
    use mpi_io_errors, only: fail_input_io_mpi
    use yes_no_input, only: getyn

    ! Procedure interfaces
    use read_field_shift_factors_interface, only: read_field_shift_factors
    use read_seltzer_field_shift_factors_interface, only: read_seltzer_field_shift_factors
    use calculate_seltzer_field_shift_factors_interface, only: calculate_seltzer_field_shift_factors
    use read_field_shift_factors_mpi_interface, only: read_field_shift_factors_mpi
    use read_seltzer_field_shift_factors_mpi_interface, only: read_seltzer_field_shift_factors_mpi
    use calculate_field_shift_factors_interface, only: calculate_field_shift_factors
    use calculate_field_shift_factors_mpi_interface, only: calculate_field_shift_factors_mpi
    use calculate_seltzer_field_shift_factors_mpi_interface, only: calculate_seltzer_field_shift_factors_mpi
    use read_specific_mass_shift_interface, only: read_specific_mass_shift
    use read_specific_mass_shift_mpi_interface, only: read_specific_mass_shift_mpi
    use calculate_specific_mass_shift_interface, only: calculate_specific_mass_shift
    use calculate_specific_mass_shift_mpi_interface, only: calculate_specific_mass_shift_mpi
    use fit_electron_density_interface, only: fit_electron_density

    ! Module dependencies
    use mpi_utilities, only: gdrsummpi_root
    use run_timer_start, only: starttime
    use run_timer_finish, only: stoptime
    implicit none (type, external)

    ! Arguments
    character(len=24), intent(in) :: name

    ! Local declarations
    real(real64), pointer :: dint1vec(:, :, :)
    real(real64), pointer   :: dens1vec(:, :)
    real(real64), pointer     :: densfit(:)
    real(real64), pointer   :: fmat(:, :)
    real(real64), pointer     :: rho(:)
    real(real64), pointer     :: res(:)
    real(real64) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
    real(real64) :: dint1(max_orbitals, max_orbitals), dint2(max_orbitals, max_orbitals), dint3(max_orbitals, max_orbitals), &
        dint4(max_orbitals, max_orbitals), dint5(max_orbitals, max_orbitals), &
        dint6(max_orbitals, max_orbitals), &
        dint7(max_orbitals, max_orbitals)
    real(real64) :: au2fm, rcre, hzsmsu, hznmsu, fo90
    logical :: avail_tb, avail_ob, yes2, saveyes, serial_save
    integer :: i, orbital_index, start_clock_count, doit_ob, doit_tb, nrnuc, yes2_mpi
    integer :: avail_ob_mpi, avail_tb_mpi, serial_save_mpi, cache_status
    if (mpi_rank == 0) then
        write(*,*) "-------------------------------"
        write(*,*) "RIS_CAL_MPI: Execution Begins ..."
        write(*,*) "-------------------------------"
    end if

    doit_ob = 0
    doit_tb = 0
    yes2_mpi = 0
    avail_ob_mpi = 0
    avail_tb_mpi = 0
    serial_save_mpi = 0

    fo90 = nuclear_model_parameters(1)+2.0_real64*log(3.0_real64)*nuclear_model_parameters(2)
    do i=1,max_radial_points
        if(radial_grid(i+1)>fo90.and.radial_grid(i)<=fo90) then
            nrnuc = i+1
            exit
        end if
    end do

    allocate( dens1vec(nvec,nrnuc),dint1vec(max_orbitals,max_orbitals,nrnuc) )
    allocate( densfit(nrnuc) )
    dint1vec(:,:,:) = 0.0_real64
    dens1vec(:,:) = 0.0_real64
    allocate( fmat(6,nvec),rho(nvec),res(nvec) )

    call alloc (smsc1, nvec,"SMSC1", "RIS_CAL_MPI")
    call alloc (smsc2, nvec,"SMSC2", "RIS_CAL_MPI")
    call alloc (dens1, nvec,"DENS1","RIS_CAL_MPI")
    call alloc (dens2, nvec,"DENS2","RIS_CAL_MPI")
    call alloc (dens3, nvec,"DENS3","RIS_CAL_MPI")
    call alloc (dens4, nvec,"DENS4","RIS_CAL_MPI")
    call alloc (dens5, nvec,"DENS5","RIS_CAL_MPI")
    call alloc (dens6, nvec,"DENS6","RIS_CAL_MPI")
    call alloc (dens7, nvec,"DENS7","RIS_CAL_MPI")

    call starttime (start_clock_count, "RIS_CAL_MPI")

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

    if (mpi_rank == 0) then
        write(*,*) " Compute higher order field shift electronic factors?"
        yes2 = getyn ()
        if (yes2) then
            yes2_mpi = 1
        else
            yes2_mpi = 0
        end if
    end if
    call mpi_bcast(yes2_mpi, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
    yes2 = yes2_mpi == 1

    do i = 1,orbital_count
        do orbital_index = 1,orbital_count
            if (orbital_kappa(i)==orbital_kappa(orbital_index)) then
                dint1(i,orbital_index) = field_shift_radial_integral(i,orbital_index)
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

    do i = 1,orbital_count
        do orbital_index = 1,orbital_count
            if (i/=orbital_index) then
                rcre = cre(orbital_kappa(i),1,orbital_kappa(orbital_index))
                if (dabs(rcre) > cutoff) then
                    vint  (i,orbital_index) = vinti_radial_integral(i,orbital_index)
                    vint2(i,orbital_index) = vint(i,orbital_index)                                  &
                        + specific_mass_shift_radial_integral_2(i,orbital_index)/rcre                                &
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

    if (mpi_rank == 0) then
        call open_isotope_shift_angular_cache(name,avail_ob,1,cache_status,one_body_cache_unit)
        if (cache_status /= 0) call fail_input_io_mpi(cache_status,one_body_cache_unit,"ris4: one-body angular cache "//trim(name))
        if (avail_ob) then
            close(one_body_cache_unit,iostat=cache_status)
            if (cache_status /= 0) call fail_input_io_mpi(cache_status,one_body_cache_unit,"ris4: close one-body angular cache")
        end if
        call open_isotope_shift_angular_cache(name,avail_tb,2,cache_status,two_body_cache_unit)
        if (cache_status /= 0) call fail_input_io_mpi(cache_status,two_body_cache_unit,"ris4: two-body angular cache "//trim(name))
        if (avail_tb) then
            close(two_body_cache_unit,iostat=cache_status)
            if (cache_status /= 0) call fail_input_io_mpi(cache_status,two_body_cache_unit,"ris4: close two-body angular cache")
        end if
        avail_ob_mpi = merge(1, 0, avail_ob)
        avail_tb_mpi = merge(1, 0, avail_tb)
        if ((.not. avail_ob) .and. (.not. avail_tb)) then
            print *," Save ang. coefficients of one- and two-body op.?"
            saveyes = getyn ()
            print *
            if(saveyes) then
                doit_ob = 1
                doit_tb = 1
            end if
        else if (.not. avail_ob) then
            print *," Save ang. coefficients of one-body op. ?"
            saveyes = getyn ()
            print *
            if(saveyes) doit_ob = 1
        else if (.not. avail_tb) then
            print *," Save ang. coefficients of two-body op. ?"
            saveyes = getyn ()
            print *
            if(saveyes) doit_tb = 1
        end if
    end if
    call mpi_bcast(avail_ob_mpi, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
    call mpi_bcast(avail_tb_mpi, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
    call mpi_bcast(doit_ob, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
    call mpi_bcast(doit_tb, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
    avail_ob = avail_ob_mpi == 1
    avail_tb = avail_tb_mpi == 1
    if (mpi_rank == 0) then
        ! Saving angular coefficients still has to stay on the root
        ! process because the serial file format is shared.
        serial_save = ((.not. avail_ob .and. doit_ob == 1) .or.      &
            (.not. avail_tb .and. doit_tb == 1))
        serial_save_mpi = merge(1, 0, serial_save)
    end if
    call mpi_bcast(serial_save_mpi, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
    serial_save = serial_save_mpi == 1

    if (serial_save) then
        call mpi_barrier(mpi_comm_world, mpi_error_code)
        if (mpi_rank == 0) then
            orbital_index = index(name," ")
            if (avail_ob) then
                open(newunit=one_body_cache_unit,file = name(1:orbital_index-1)//".IOB",status="OLD",form="UNFORMATTED")
                if (yes2) then
                    call read_seltzer_field_shift_factors(dint1,dint2,dint3,dint4,dint5,dint6, &
                        dint7,dint1vec,dens1vec,nrnuc)
                else
                    call read_field_shift_factors(dint1,dint2,dint3,dint4,dint5,dint6,dint7)
                end if
                close(one_body_cache_unit)
            else
                if (doit_ob == 1) then
                    open(newunit=one_body_cache_unit,file = name(1:orbital_index-1)//".IOB",status="UNKNOWN",form="UNFORMATTED")
                end if
                if (yes2) then
                    call calculate_seltzer_field_shift_factors(doit_ob,dint1,dint2,dint3,dint4,dint5,dint6, &
                        dint7,dint1vec,dens1vec,nrnuc)
                else
                    call calculate_field_shift_factors(doit_ob,dint1,dint2,dint3,dint4,dint5,dint6,dint7)
                end if
            end if

            if (avail_tb) then
                open(newunit=two_body_cache_unit,file = name(1:orbital_index-1)//".ITB",status="OLD",form="UNFORMATTED")
                call read_specific_mass_shift(vint,vint2)
                close(two_body_cache_unit)
            else
                if (doit_tb == 1) then
                    open(newunit=two_body_cache_unit,file = name(1:orbital_index-1)//".ITB",status="UNKNOWN",form="UNFORMATTED")
                end if
                call calculate_specific_mass_shift(doit_tb,vint,vint2)
            end if
        end if
        call mpi_barrier(mpi_comm_world, mpi_error_code)
    else
        if (avail_ob) then
            orbital_index = index(name," ")
            open(newunit=one_body_cache_unit,file = name(1:orbital_index-1)//".IOB",status="OLD",form="UNFORMATTED")
            if (yes2) then
                call read_seltzer_field_shift_factors_mpi(dint1,dint2,dint3,dint4,dint5,dint6, &
                    dint7,dint1vec,dens1vec,nrnuc)
            else
                call read_field_shift_factors_mpi(dint1,dint2,dint3,dint4,dint5,dint6,dint7)
            end if
            close(one_body_cache_unit)
        else
            if(yes2) then
                call calculate_seltzer_field_shift_factors_mpi(doit_ob,dint1,dint2,dint3,dint4,dint5,dint6, &
                    dint7,dint1vec,dens1vec,nrnuc)
            else
                call calculate_field_shift_factors_mpi(doit_ob,dint1,dint2,dint3,dint4,dint5,dint6,dint7)
            end if
        end if

        if (avail_tb) then
            orbital_index = index(name," ")
            open(newunit=two_body_cache_unit,file = name(1:orbital_index-1)//".ITB",status="OLD",form="UNFORMATTED")
            call read_specific_mass_shift_mpi(vint,vint2)
            close(two_body_cache_unit)
        else
            call calculate_specific_mass_shift_mpi(doit_tb,vint,vint2)
        end if

        call mpi_barrier(mpi_comm_world, mpi_error_code)

        call gdrsummpi_root(smsc1(:), nvec)
        call gdrsummpi_root(smsc2(:), nvec)
        call gdrsummpi_root(dens1(:), nvec)
        call gdrsummpi_root(dens2(:), nvec)
        call gdrsummpi_root(dens3(:), nvec)
        call gdrsummpi_root(dens4(:), nvec)
        call gdrsummpi_root(dens5(:), nvec)
        call gdrsummpi_root(dens6(:), nvec)
        call gdrsummpi_root(dens7(:), nvec)
        if (yes2) call gdrsummpi_root(dens1vec(:,:), nvec*nrnuc)
    end if

    if (mpi_rank == 0) then
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
            write (summary_output_unit,'(20x,3d20.10,2x,"(GHz u)")') hznmsu*(1.0e-9_real64), &
                (dens2(i)-dens7(i))*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64), &
                dens2(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64)
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
                call fit_electron_density(radial_grid,dens1vec(i,:),nuclear_charge,nuclear_model_parameters,nrnuc,fmat(:,i), &
                    & rho(i), &
                    & res(i))
                write (summary_output_unit,'(1x,i3,5x,2a4,3x,d20.10)') ivec(i),right_justified_j_labels(iatjpo(i)), &
                    parity_labels((iaspar(i)+3)/2),rho(i)
            end do
            write (summary_output_unit, &
                & '(//" Level  J Parity  Field shift electronic factors and average point discrepancy in fit"//&
                &24x,"F0 (GHz/fm^2)",6x," F2 (GHz/fm^4)",6x," F4 (GHz/fm^6)",6x," F6 (GHz/fm^8)",6x,&
                &" Disc. (per mille)")')
            do i = 1,nvec
                write (summary_output_unit,'(1x,i3,5x,2a4,3x,4d20.10,f10.4)') ivec(i),right_justified_j_labels(iatjpo(i)), &
                    parity_labels((iaspar(i)+3)/2),fmat(1,i),fmat(2,i), &
                    fmat(3,i),fmat(4,i),res(i)
            end do
            write (summary_output_unit,'(//&
                &" Level  J Parity  Field shift electronic factors (corrected for varying density inside nucleus)"//&
                &24x,"F0VED0 (GHz/fm^2)   F0VED1 (GHz/fm^4)")')
            do i = 1,nvec
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
    end if

    call dalloc (smsc1,"SMSC1","RIS_CAL_MPI")
    call dalloc (smsc2,"SMSC2","RIS_CAL_MPI")
    call dalloc (dens1,"DENS1","RIS_CAL_MPI")
    call dalloc (dens2,"DENS2","RIS_CAL_MPI")
    call dalloc (dens3,"DENS3","RIS_CAL_MPI")
    call dalloc (dens4,"DENS4","RIS_CAL_MPI")
    call dalloc (dens5,"DENS5","RIS_CAL_MPI")
    call dalloc (dens6,"DENS6","RIS_CAL_MPI")
    call dalloc (dens7,"DENS7","RIS_CAL_MPI")

    if (mpi_rank == 0) then
        write(*,*) "-------------------------------"
        write(*,*) "RIS_CAL_MPI: Execution Finished ..."
        write(*,*) "-------------------------------"
    end if
    call stoptime (start_clock_count, "RIS_CAL_MPI")
    return
    return
end subroutine calculate_isotope_shifts_mpi
