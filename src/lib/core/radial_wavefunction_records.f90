! Match radial records to the orbital basis by n and kappa, then interpolate and orthogonalize the loaded
! orbitals.
module radial_wavefunction_records
    implicit none (type, external)

    ! Module interface
    private
    public :: read_radial_wavefunction_records
contains
    subroutine read_radial_wavefunction_records(wavefunction_input_unit, load_status)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_radial_points
        use array_memory_management, only: alloc, dalloc
        use debug_output_settings, only: rmcdhf_debug_flags
        use atomic_calculation_state,           only: speed_of_light_au, nuclear_charge
        use radial_grid_state, only: radial_point_count
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count
        use orbital_csf_basis, only: orbital_count, orbital_energy_parameters, orbital_kappa, orbital_origin_exponents, &
            orbital_principal_n, orbital_symmetry_labels
        use radial_wavefunctions,          only: orbital_origin_parameters, radial_large_components, radial_small_components
        use radial_orbital_interpolation, only: interpolate_radial_orbital
        use radial_orbital_orthogonalization, only: orthogonalize_radial_orbitals
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: wavefunction_input_unit
        integer, intent(out) :: load_status

        ! Local declarations
        integer :: orbital_index, radial_index, abs_orbital_kappa, loaded_orbital_count, read_status, file_principal_n, &
            & file_kappa, file_grid_size
        real(real64) :: charge_times_alpha_squared, kappa_squared, file_energy_parameter, file_origin_coefficient, &
            & original_norm_squared
        real(real64), pointer :: file_large_component(:), file_small_component(:), file_radial_grid(:)
        write (output_unit, *) "Loading Radial WaveFunction File ..."
        call alloc(radial_large_components, max_radial_points, orbital_count, "radial_large_components", &
            "read_radial_wavefunction_records")
        call alloc(radial_small_components, max_radial_points, orbital_count, "radial_small_components", &
            "read_radial_wavefunction_records")

        charge_times_alpha_squared = nuclear_charge/speed_of_light_au
        charge_times_alpha_squared = charge_times_alpha_squared*charge_times_alpha_squared
        do orbital_index = 1, orbital_count
            radial_large_components(:radial_point_count,orbital_index) = 0.0_real64
            radial_small_components(:radial_point_count,orbital_index) = 0.0_real64
            orbital_energy_parameters(orbital_index) = -1.0_real64
            abs_orbital_kappa = abs(orbital_kappa(orbital_index))
            if (nuclear_model_parameter_count > 0) then
                orbital_origin_exponents(orbital_index) = dble(abs_orbital_kappa)
            else if (nuclear_model_parameter_count == 0) then
                kappa_squared = dble(abs_orbital_kappa*abs_orbital_kappa)
                if (kappa_squared >= charge_times_alpha_squared) then
                    orbital_origin_exponents(orbital_index) = sqrt(kappa_squared - charge_times_alpha_squared)
                else
                    stop "read_radial_wavefunction_records: Inappropriate gamma"
                end if
            end if
        end do
        if (rmcdhf_debug_flags(3)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,"From SUBROUTINE LODRWF:"/," Orbital",8x,"Eigenvalue",19x,"Norm")')
        end if
        loaded_orbital_count = 0
        ! Load each binary record in file order and install it only in matching orbital slots.
        do
            read (wavefunction_input_unit, iostat=read_status) file_principal_n, file_kappa, file_energy_parameter, file_grid_size
            if (read_status /= 0) exit
            call alloc (file_large_component,file_grid_size, "file_large_component", "read_radial_wavefunction_records")
            call alloc (file_small_component,file_grid_size, "file_small_component", "read_radial_wavefunction_records")
            call alloc (file_radial_grid,file_grid_size, "file_radial_grid", "read_radial_wavefunction_records")

            ! The binary format stores the coefficient and both components together, then radii separately.
            read (wavefunction_input_unit) file_origin_coefficient, (file_large_component(radial_index),radial_index=1, &
                file_grid_size), &
                & (file_small_component(radial_index),radial_index=1,file_grid_size)
            read (wavefunction_input_unit) (file_radial_grid(radial_index),radial_index=1,file_grid_size)

            do orbital_index = 1, orbital_count
                if (.not.(orbital_energy_parameters(orbital_index)<0.0_real64 .and. &
                    & file_principal_n==orbital_principal_n(orbital_index) .and. file_kappa==orbital_kappa(orbital_index))) cycle
                orbital_origin_parameters(orbital_index) = file_origin_coefficient
                orbital_energy_parameters(orbital_index) = file_energy_parameter
                call interpolate_radial_orbital(file_large_component, file_small_component, file_grid_size, file_radial_grid, &
                    & orbital_index, &
                    & original_norm_squared)
                if (rmcdhf_debug_flags(3)) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)') &
                        orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), &
                        & orbital_energy_parameters(orbital_index), original_norm_squared
                end if
                loaded_orbital_count = loaded_orbital_count + 1
            end do
            call dalloc (file_large_component, "file_large_component", "read_radial_wavefunction_records" )
            call dalloc (file_small_component, "file_small_component", "read_radial_wavefunction_records")
            call dalloc (file_radial_grid, "file_radial_grid", "read_radial_wavefunction_records")
        end do
        if (rmcdhf_debug_flags(3)) then
            call ensure_debug_output()
            write (angular_debug_unit, *) " orbitals renormalised;"
        end if
        if (loaded_orbital_count < orbital_count) then
            load_status = 1
            return
        end if
        call orthogonalize_radial_orbitals
        if (rmcdhf_debug_flags(3)) then
            call ensure_debug_output()
            write (angular_debug_unit, *) " orbitals orthogonalised and renormalised;"
        end if
        load_status = 0
    end subroutine read_radial_wavefunction_records
end module radial_wavefunction_records
