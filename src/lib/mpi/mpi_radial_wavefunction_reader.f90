module mpi_radial_wavefunction_reader
    implicit none (type, external)

    ! Module interface
    private
    public :: read_mpi_radial_wavefunction_records
contains
    ! Read orbital records on rank 0 and interpolate the broadcast wavefunctions on each rank.
    subroutine read_mpi_radial_wavefunction_records (wavefunction_input_unit, load_status)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output

        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only:  max_radial_points
        use array_memory_management, only: alloc, dalloc
        use debug_output_settings, only: rmcdhf_debug_flags
        use atomic_calculation_state, only: speed_of_light_au, nuclear_charge
        use radial_grid_state, only: radial_point_count
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count
        use orbital_csf_basis, only: orbital_count, orbital_energy_parameters, orbital_kappa, orbital_origin_exponents, &
            orbital_principal_n, orbital_symmetry_labels
        use radial_wavefunctions, only: orbital_origin_parameters, radial_large_components, radial_small_components
        use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_integer, mpi_rank
        use radial_orbital_interpolation, only: interpolate_radial_orbital
        use radial_orbital_orthogonalization, only: orthogonalize_radial_orbitals
        use mpi_utilities, only: finalize_mpi_runtime
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: wavefunction_input_unit
        integer, intent(out) :: load_status

        ! Local declarations
        integer :: orbital_index, radial_index, abs_orbital_kappa, loaded_orbital_count
        integer :: read_status, file_principal_n, file_kappa, file_grid_size
        real(real64) :: squared_charge_to_light_speed_ratio, kappa_squared
        real(real64) :: file_energy_parameter, file_origin_coefficient, original_norm_squared
        real(real64), pointer :: file_large_component(:), file_small_component(:), file_radial_grid(:)
        call alloc(radial_large_components, max_radial_points, orbital_count, "PF", "LODRWFMPI")
        call alloc(radial_small_components, max_radial_points, orbital_count, "QF", "LODRWFMPI")

        ! Prepare the point-nucleus origin exponent. Finite nuclei use the absolute value of kappa.
        squared_charge_to_light_speed_ratio = nuclear_charge / speed_of_light_au
        squared_charge_to_light_speed_ratio = squared_charge_to_light_speed_ratio * squared_charge_to_light_speed_ratio

        do orbital_index = 1, orbital_count

            do radial_index = 1, radial_point_count
                radial_large_components(radial_index,orbital_index) = 0.0_real64
                radial_small_components(radial_index,orbital_index) = 0.0_real64
            end do

            ! A negative energy parameter marks an orbit still eligible for matching.
            orbital_energy_parameters(orbital_index) = -1.0_real64
            abs_orbital_kappa = abs (orbital_kappa(orbital_index))
            if (nuclear_model_parameter_count > 0) then
                orbital_origin_exponents(orbital_index) = dble (abs_orbital_kappa)
            else if (nuclear_model_parameter_count == 0) then
                kappa_squared = dble (abs_orbital_kappa*abs_orbital_kappa)
                if (kappa_squared >= squared_charge_to_light_speed_ratio) then
                    orbital_origin_exponents(orbital_index) = sqrt (kappa_squared-squared_charge_to_light_speed_ratio)
                else
                    call finalize_mpi_runtime ("lodrwfmpi: Inappropriate gamma", mpi_rank)
                end if
            end if
        end do
        if (rmcdhf_debug_flags(3) .and. mpi_rank == 0) then
            call ensure_debug_output()
            write (angular_debug_unit,'(/"From SUBROUTINE LODRWFmpi:"/" Orbital",8x,&
                &"Eigenvalue",19x,"Norm")')
        end if
        loaded_orbital_count = 0
        ! Each orbital starts with a record containing its identifier, energy parameter and grid size.
        do
            if (mpi_rank == 0) read (wavefunction_input_unit,iostat = read_status) file_principal_n,file_kappa, &
                file_energy_parameter,file_grid_size
            ! Share the read status so all ranks take the same branch before allocating and receiving the payload.
            call mpi_bcast (read_status,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
            call mpi_bcast (file_principal_n,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
            call mpi_bcast (file_kappa,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
            call mpi_bcast (file_grid_size,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
            call mpi_bcast (file_energy_parameter,1,mpi_double_precision,0,mpi_comm_world,mpi_error_code)
            if (read_status /= 0) exit

            call alloc (file_large_component,file_grid_size, "PA", "LODRWFMPI")
            call alloc (file_small_component,file_grid_size, "QA", "LODRWFMPI")
            call alloc (file_radial_grid,file_grid_size, "RA", "LODRWFMPI")

            if (mpi_rank == 0) then
                read (wavefunction_input_unit) file_origin_coefficient, (file_large_component(radial_index), radial_index = 1, &
                    file_grid_size), &
                    (file_small_component(radial_index), radial_index = 1,file_grid_size)
                read (wavefunction_input_unit) (file_radial_grid(radial_index),radial_index = 1,file_grid_size)
            end if

            call mpi_bcast (file_origin_coefficient,1,mpi_double_precision,0, &
                mpi_comm_world,mpi_error_code)
            call mpi_bcast (file_large_component,file_grid_size,mpi_double_precision,0, &
                mpi_comm_world,mpi_error_code)
            call mpi_bcast (file_small_component,file_grid_size,mpi_double_precision,0, &
                mpi_comm_world,mpi_error_code)
            call mpi_bcast (file_radial_grid,file_grid_size,mpi_double_precision,0, &
                mpi_comm_world,mpi_error_code)

            do orbital_index = 1, orbital_count
                if (.not.(orbital_energy_parameters(orbital_index) < 0.0_real64 .and. &
                    file_principal_n == orbital_principal_n(orbital_index) .and. &
                    file_kappa == orbital_kappa(orbital_index))) cycle
                orbital_origin_parameters(orbital_index) = file_origin_coefficient
                orbital_energy_parameters(orbital_index) = file_energy_parameter
                call interpolate_radial_orbital (file_large_component, file_small_component, file_grid_size, file_radial_grid, &
                    orbital_index, original_norm_squared)
                if (rmcdhf_debug_flags(3) .and. mpi_rank == 0) then
                    call ensure_debug_output()
                    write (angular_debug_unit,'(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)') &
                        orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), &
                        orbital_energy_parameters(orbital_index), original_norm_squared
                end if
                loaded_orbital_count = loaded_orbital_count + 1
            end do

            call dalloc (file_large_component, "PA", "LODRWFMPI")
            call dalloc (file_small_component, "QA", "LODRWFMPI")
            call dalloc (file_radial_grid, "RA", "LODRWFMPI")
        end do
        if (rmcdhf_debug_flags(3) .and. mpi_rank == 0) then
            call ensure_debug_output()
            write (angular_debug_unit,*) " orbitals renormalised;"
        end if
        ! Return failure if fewer orbital matches were loaded than requested.
        if (loaded_orbital_count < orbital_count) then
            load_status = 1
            return
        end if
        ! Orthogonalize loaded orbitals within each kappa symmetry.
        call orthogonalize_radial_orbitals

        if (rmcdhf_debug_flags(3) .and. mpi_rank == 0) then
            call ensure_debug_output()
            write (angular_debug_unit,*) "orbitals orthogonalised and renormalised;"
        end if
        load_status = 0

        return
    end subroutine read_mpi_radial_wavefunction_records
end module mpi_radial_wavefunction_reader
