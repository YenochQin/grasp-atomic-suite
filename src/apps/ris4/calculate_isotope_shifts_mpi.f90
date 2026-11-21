module calculate_isotope_shifts_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_isotope_shifts_mpi
contains

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
        use debug_output_settings, only: cutoff
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
        use angular_symmetry_labels, only: right_justified_j_labels, parity_labels
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues
        use mpi_state, only: mpi_barrier, mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_rank

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
        use mpi_io_errors, only: fail_input_io_mpi
        use yes_no_input, only: read_yes_no

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
        use mpi_utilities, only: reduce_mpi_real_sum_to_root
        use run_timer_start, only: begin_run_timer
        use run_timer_finish, only: finish_run_timer
        implicit none (type, external)

        ! Arguments
        character(len=24), intent(in) :: name

        ! Local declarations
        integer :: first_name_blank
        real(real64), pointer :: field_shift_radial_integral_profiles(:, :, :)
        real(real64), pointer   :: level_electron_density_profiles(:, :)
        real(real64), pointer     :: electron_density_fit(:)
        real(real64), pointer   :: density_fit_parameters(:, :)
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
        logical :: avail_tb, avail_ob, yes2, saveyes, serial_save
        integer :: i, orbital_index, start_clock_count, one_body_cache_mode, two_body_cache_mode, &
            & nuclear_density_grid_point_count, yes2_mpi
        integer :: avail_ob_mpi, avail_tb_mpi, serial_save_mpi, cache_status
        if (mpi_rank == 0) then
            write(*,*) "-------------------------------"
            write(*,*) "calculate_isotope_shifts_mpi: Execution Begins ..."
            write(*,*) "-------------------------------"
        end if

        one_body_cache_mode = 0
        two_body_cache_mode = 0
        yes2_mpi = 0
        avail_ob_mpi = 0
        avail_tb_mpi = 0
        serial_save_mpi = 0

        nuclear_density_fit_radius = nuclear_model_parameters(1)+2.0_real64*log(3.0_real64)*nuclear_model_parameters(2)
        do i=1,max_radial_points
            if(radial_grid(i+1)>nuclear_density_fit_radius.and.radial_grid(i)<=nuclear_density_fit_radius) then
                nuclear_density_grid_point_count = i+1
                exit
            end if
        end do

        allocate( level_electron_density_profiles(active_level_count,nuclear_density_grid_point_count), &
            & field_shift_radial_integral_profiles(max_orbitals,max_orbitals,nuclear_density_grid_point_count) )
        allocate( electron_density_fit(nuclear_density_grid_point_count) )
        field_shift_radial_integral_profiles(:,:,:) = 0.0_real64
        level_electron_density_profiles(:,:) = 0.0_real64
        allocate( density_fit_parameters(6,active_level_count),effective_electron_density(active_level_count), &
            density_fit_residual(active_level_count) )

        call alloc (leading_sms_by_level, active_level_count,"leading_sms_by_level", "calculate_isotope_shifts_mpi")
        call alloc (total_sms_by_level, active_level_count,"total_sms_by_level", "calculate_isotope_shifts_mpi")
        call alloc (contact_density_by_level, active_level_count,"contact_density_by_level","calculate_isotope_shifts_mpi")
        call alloc (total_nms_by_level, active_level_count,"total_nms_by_level","calculate_isotope_shifts_mpi")
        call alloc (radius_by_level, active_level_count,"radius_by_level","calculate_isotope_shifts_mpi")
        call alloc (radius_squared_by_level, active_level_count,"radius_squared_by_level","calculate_isotope_shifts_mpi")
        call alloc (inverse_radius_by_level, active_level_count,"inverse_radius_by_level","calculate_isotope_shifts_mpi")
        call alloc (inverse_radius_squared_by_level, active_level_count,"inverse_radius_squared_by_level", &
            "calculate_isotope_shifts_mpi")
        call alloc (leading_nms_by_level, active_level_count,"leading_nms_by_level","calculate_isotope_shifts_mpi")

        call begin_run_timer (start_clock_count, "RIS_CAL_MPI")

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

        if (mpi_rank == 0) then
            write(*,*) " Compute higher order field shift electronic factors?"
            yes2 = read_yes_no ()
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
                    contact_density_integrals(i,orbital_index) = field_shift_radial_integral(i,orbital_index)
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

        do i = 1,orbital_count
            do orbital_index = 1,orbital_count
                if (i/=orbital_index) then
                    reduced_rank1_angular_element = reduced_angular_element(orbital_kappa(i),1,orbital_kappa(orbital_index))
                    if (dabs(reduced_rank1_angular_element) > cutoff) then
                        vinti_integrals  (i,orbital_index) = vinti_radial_integral(i,orbital_index)
                        recoil_corrected_vinti_integrals(i,orbital_index) = vinti_integrals(i, &
                            & orbital_index)                                  &
                            + specific_mass_shift_radial_integral_2(i,orbital_index)/ &
                                & reduced_rank1_angular_element                                &
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

        if (mpi_rank == 0) then
            call open_isotope_shift_angular_cache(name,avail_ob,1,cache_status,one_body_cache_unit)
            if (cache_status /= 0) call fail_input_io_mpi(cache_status,one_body_cache_unit, &
                & "ris4: one-body angular cache "//trim(name))
            if (avail_ob) then
                close(one_body_cache_unit,iostat=cache_status)
                if (cache_status /= 0) call fail_input_io_mpi(cache_status,one_body_cache_unit,"ris4: close one-body angular cache")
            end if
            call open_isotope_shift_angular_cache(name,avail_tb,2,cache_status,two_body_cache_unit)
            if (cache_status /= 0) call fail_input_io_mpi(cache_status,two_body_cache_unit, &
                & "ris4: two-body angular cache "//trim(name))
            if (avail_tb) then
                close(two_body_cache_unit,iostat=cache_status)
                if (cache_status /= 0) call fail_input_io_mpi(cache_status,two_body_cache_unit,"ris4: close two-body angular cache")
            end if
            avail_ob_mpi = merge(1, 0, avail_ob)
            avail_tb_mpi = merge(1, 0, avail_tb)
            if ((.not. avail_ob) .and. (.not. avail_tb)) then
                print *," Save ang. coefficients of one- and two-body op.?"
                saveyes = read_yes_no ()
                print *
                if(saveyes) then
                    one_body_cache_mode = 1
                    two_body_cache_mode = 1
                end if
            else if (.not. avail_ob) then
                print *," Save ang. coefficients of one-body op. ?"
                saveyes = read_yes_no ()
                print *
                if(saveyes) one_body_cache_mode = 1
            else if (.not. avail_tb) then
                print *," Save ang. coefficients of two-body op. ?"
                saveyes = read_yes_no ()
                print *
                if(saveyes) two_body_cache_mode = 1
            end if
        end if
        call mpi_bcast(avail_ob_mpi, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
        call mpi_bcast(avail_tb_mpi, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
        call mpi_bcast(one_body_cache_mode, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
        call mpi_bcast(two_body_cache_mode, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
        avail_ob = avail_ob_mpi == 1
        avail_tb = avail_tb_mpi == 1
        if (mpi_rank == 0) then
            ! Saving angular coefficients still has to stay on the root
            ! process because the serial file format is shared.
            serial_save = ((.not. avail_ob .and. one_body_cache_mode == 1) .or.      &
                (.not. avail_tb .and. two_body_cache_mode == 1))
            serial_save_mpi = merge(1, 0, serial_save)
        end if
        call mpi_bcast(serial_save_mpi, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
        serial_save = serial_save_mpi == 1

        if (serial_save) then
            call mpi_barrier(mpi_comm_world, mpi_error_code)
            if (mpi_rank == 0) then
                first_name_blank = index(name," ")
                if (avail_ob) then
                    open(newunit=one_body_cache_unit,file = name(1:first_name_blank-1)//".IOB",status="OLD",form="UNFORMATTED")
                    if (yes2) then
                        call read_seltzer_field_shift_factors(contact_density_integrals,total_nms_integrals, &
                            & radius_integrals,radius_squared_integrals,inverse_radius_integrals,inverse_radius_squared_integrals, &
                            leading_nms_integrals,field_shift_radial_integral_profiles,level_electron_density_profiles, &
                                & nuclear_density_grid_point_count)
                    else
                        call read_field_shift_factors(contact_density_integrals,total_nms_integrals,radius_integrals, &
                            & radius_squared_integrals,inverse_radius_integrals,inverse_radius_squared_integrals, &
                            & leading_nms_integrals)
                    end if
                    close(one_body_cache_unit)
                else
                    if (one_body_cache_mode == 1) then
                        open(newunit=one_body_cache_unit,file = name(1:first_name_blank-1)//".IOB",status="UNKNOWN",form= &
                            & "UNFORMATTED")
                    end if
                    if (yes2) then
                        call calculate_seltzer_field_shift_factors(one_body_cache_mode,contact_density_integrals, &
                            & total_nms_integrals,radius_integrals,radius_squared_integrals,inverse_radius_integrals, &
                            & inverse_radius_squared_integrals, &
                            leading_nms_integrals,field_shift_radial_integral_profiles,level_electron_density_profiles, &
                                & nuclear_density_grid_point_count)
                    else
                        call calculate_field_shift_factors(one_body_cache_mode,contact_density_integrals, &
                            & total_nms_integrals,radius_integrals,radius_squared_integrals,inverse_radius_integrals, &
                            & inverse_radius_squared_integrals,leading_nms_integrals)
                    end if
                end if

                if (avail_tb) then
                    open(newunit=two_body_cache_unit,file = name(1:first_name_blank-1)//".ITB",status="OLD",form="UNFORMATTED")
                    call read_specific_mass_shift(vinti_integrals,recoil_corrected_vinti_integrals)
                    close(two_body_cache_unit)
                else
                    if (two_body_cache_mode == 1) then
                        open(newunit=two_body_cache_unit,file = name(1:first_name_blank-1)//".ITB",status="UNKNOWN",form= &
                            & "UNFORMATTED")
                    end if
                    call calculate_specific_mass_shift(two_body_cache_mode,vinti_integrals,recoil_corrected_vinti_integrals)
                end if
            end if
            call mpi_barrier(mpi_comm_world, mpi_error_code)
        else
            if (avail_ob) then
                first_name_blank = index(name," ")
                open(newunit=one_body_cache_unit,file = name(1:first_name_blank-1)//".IOB",status="OLD",form="UNFORMATTED")
                if (yes2) then
                    call read_seltzer_field_shift_factors_mpi(contact_density_integrals,total_nms_integrals, &
                        & radius_integrals,radius_squared_integrals,inverse_radius_integrals,inverse_radius_squared_integrals, &
                        leading_nms_integrals,field_shift_radial_integral_profiles,level_electron_density_profiles, &
                            & nuclear_density_grid_point_count)
                else
                    call read_field_shift_factors_mpi(contact_density_integrals,total_nms_integrals,radius_integrals, &
                        & radius_squared_integrals,inverse_radius_integrals,inverse_radius_squared_integrals,leading_nms_integrals)
                end if
                close(one_body_cache_unit)
            else
                if(yes2) then
                    call calculate_seltzer_field_shift_factors_mpi(one_body_cache_mode,contact_density_integrals, &
                        & total_nms_integrals,radius_integrals,radius_squared_integrals,inverse_radius_integrals, &
                        & inverse_radius_squared_integrals, &
                        leading_nms_integrals,field_shift_radial_integral_profiles,level_electron_density_profiles, &
                            & nuclear_density_grid_point_count)
                else
                    call calculate_field_shift_factors_mpi(one_body_cache_mode,contact_density_integrals, &
                        & total_nms_integrals,radius_integrals,radius_squared_integrals,inverse_radius_integrals, &
                        & inverse_radius_squared_integrals,leading_nms_integrals)
                end if
            end if

            if (avail_tb) then
                first_name_blank = index(name," ")
                open(newunit=two_body_cache_unit,file = name(1:first_name_blank-1)//".ITB",status="OLD",form="UNFORMATTED")
                call read_specific_mass_shift_mpi(vinti_integrals,recoil_corrected_vinti_integrals)
                close(two_body_cache_unit)
            else
                call calculate_specific_mass_shift_mpi(two_body_cache_mode,vinti_integrals,recoil_corrected_vinti_integrals)
            end if

            call mpi_barrier(mpi_comm_world, mpi_error_code)

            call reduce_mpi_real_sum_to_root(leading_sms_by_level(:), active_level_count)
            call reduce_mpi_real_sum_to_root(total_sms_by_level(:), active_level_count)
            call reduce_mpi_real_sum_to_root(contact_density_by_level(:), active_level_count)
            call reduce_mpi_real_sum_to_root(total_nms_by_level(:), active_level_count)
            call reduce_mpi_real_sum_to_root(radius_by_level(:), active_level_count)
            call reduce_mpi_real_sum_to_root(radius_squared_by_level(:), active_level_count)
            call reduce_mpi_real_sum_to_root(inverse_radius_by_level(:), active_level_count)
            call reduce_mpi_real_sum_to_root(inverse_radius_squared_by_level(:), active_level_count)
            call reduce_mpi_real_sum_to_root(leading_nms_by_level(:), active_level_count)
            if (yes2) call reduce_mpi_real_sum_to_root(level_electron_density_profiles(:,:), active_level_count* &
                nuclear_density_grid_point_count)
        end if

        if (mpi_rank == 0) then
            write (summary_output_unit,"(a24,i3)") "Number of eigenvalues: ",active_level_count
            write(summary_output_unit,*)
            write(summary_output_unit,*)
            write (summary_output_unit,'(" Level  J Parity  Energy")')
            do i = 1, active_level_count
                write (summary_output_unit,'(1x,i3,5x,2a4,3x,d20.10,"  (a.u.)")') active_level_indices(i), &
                    right_justified_j_labels(level_two_j_plus_one(i)), &
                    & parity_labels((level_parity_signs(i)+3)/2), &
                    hamiltonian_energy_shift+shifted_eigenvalues(i)
            end do
            write (summary_output_unit,'(//" Level  J Parity  Normal mass shift parameter")')
            do i = 1,active_level_count
                write(summary_output_unit,'(/,29x,"<K^1>",13x,"<K^2+K^3>",9x,"<K^1+K^2+K^3>")')
                write (summary_output_unit,'(1x,i3,5x,2a4,3x,3d20.10,"  (a.u.)")') active_level_indices(i), &
                    right_justified_j_labels(level_two_j_plus_one(i)), &
                    & parity_labels((level_parity_signs(i)+3)/2), &
                    leading_nms_by_level(i),(total_nms_by_level(i)-leading_nms_by_level(i)),total_nms_by_level(i)
                normal_mass_shift_hz_u=leading_nms_by_level(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s
                write (summary_output_unit,'(20x,3d20.10,2x,"(GHz u)")') normal_mass_shift_hz_u*(1.0e-9_real64), &
                    (total_nms_by_level(i)-leading_nms_by_level(i))*electron_mass_u*hartree_to_wavenumber_cm* &
                        & speed_of_light_cm_per_s*(1.0e-9_real64), &
                    total_nms_by_level(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64)
            end do

            write (summary_output_unit,'(//" Level  J Parity  Specific mass shift parameter")')
            do i = 1,active_level_count
                write(summary_output_unit,'(/,29x,"<K^1>",13x,"<K^2+K^3>",9x,"<K^1+K^2+K^3>")')
                write (summary_output_unit,'(1x,i3,5x,2a4,3x,3d20.10,"  (a.u.)")') active_level_indices(i), &
                    right_justified_j_labels(level_two_j_plus_one(i)), &
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
                        & density_fit_residual(i))
                    write (summary_output_unit,'(1x,i3,5x,2a4,3x,d20.10)') active_level_indices(i), &
                        right_justified_j_labels(level_two_j_plus_one(i)), &
                        parity_labels((level_parity_signs(i)+3)/2),effective_electron_density(i)
                end do
                write (summary_output_unit, &
                    & '(//" Level  J Parity  Field shift electronic factors and average point discrepancy in fit"//&
                    &24x,"F0 (GHz/fm^2)",6x," F2 (GHz/fm^4)",6x," F4 (GHz/fm^6)",6x," F6 (GHz/fm^8)",6x,&
                    &" Disc. (per mille)")')
                do i = 1,active_level_count
                    write (summary_output_unit,'(1x,i3,5x,2a4,3x,4d20.10,f10.4)') active_level_indices(i), &
                        right_justified_j_labels(level_two_j_plus_one(i)), &
                        parity_labels((level_parity_signs(i)+3)/2),density_fit_parameters(1,i),density_fit_parameters(2,i), &
                        density_fit_parameters(3,i),density_fit_parameters(4,i),density_fit_residual(i)
                end do
                write (summary_output_unit,'(//&
                    &" Level  J Parity  Field shift electronic factors (corrected for varying density inside nucleus)"//&
                    &24x,"F0VED0 (GHz/fm^2)   F0VED1 (GHz/fm^4)")')
                do i = 1,active_level_count
                    write (summary_output_unit,'(1x,i3,5x,2a4,3x,2d20.10)') active_level_indices(i), &
                        right_justified_j_labels(level_two_j_plus_one(i)), &
                        parity_labels((level_parity_signs(i)+3)/2),density_fit_parameters(5,i),density_fit_parameters(6,i)
                end do
                write(summary_output_unit,*)
            else
                write (summary_output_unit,'(//" Electron density in atomic units"//" Level  J Parity",8x,"DENS (a.u.)"/)')
                do i = 1,active_level_count
                    write (summary_output_unit,'(1x,i3,5x,2a4,3x,d20.10)') active_level_indices(i), &
                        right_justified_j_labels(level_two_j_plus_one(i)), &
                        parity_labels((level_parity_signs(i)+3)/2),contact_density_by_level(i)
                end do
            end if
        end if

        call dalloc (leading_sms_by_level,"leading_sms_by_level","calculate_isotope_shifts_mpi")
        call dalloc (total_sms_by_level,"total_sms_by_level","calculate_isotope_shifts_mpi")
        call dalloc (contact_density_by_level,"contact_density_by_level","calculate_isotope_shifts_mpi")
        call dalloc (total_nms_by_level,"total_nms_by_level","calculate_isotope_shifts_mpi")
        call dalloc (radius_by_level,"radius_by_level","calculate_isotope_shifts_mpi")
        call dalloc (radius_squared_by_level,"radius_squared_by_level","calculate_isotope_shifts_mpi")
        call dalloc (inverse_radius_by_level,"inverse_radius_by_level","calculate_isotope_shifts_mpi")
        call dalloc (inverse_radius_squared_by_level,"inverse_radius_squared_by_level","calculate_isotope_shifts_mpi")
        call dalloc (leading_nms_by_level,"leading_nms_by_level","calculate_isotope_shifts_mpi")

        if (mpi_rank == 0) then
            write(*,*) "-------------------------------"
            write(*,*) "calculate_isotope_shifts_mpi: Execution Finished ..."
            write(*,*) "-------------------------------"
        end if
        call finish_run_timer (start_clock_count, "RIS_CAL_MPI")
        return
        return
    end subroutine calculate_isotope_shifts_mpi
end module calculate_isotope_shifts_mpi_interface
