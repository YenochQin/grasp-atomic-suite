! Match radial records to the orbital basis by n and kappa, then interpolate and orthogonalize the loaded
! orbitals.
module radial_wavefunction_records
    implicit none (type, external)

    ! Module interface
    private
    public :: lodrwf
contains
    subroutine lodrwf(wavefunction_input_unit, load_status)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_radial_points
        use array_memory_management, only: alloc, dalloc
        use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
        use atomic_calculation_state,           only: speed_of_light_au, nuclear_charge
        use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
            radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
        use orbital_csf_basis, only: orbital_count, orbital_energy_parameters, orbital_kappa, orbital_origin_exponents, &
            orbital_principal_n, orbital_symmetry_labels
        use radial_wavefunctions,          only: orbital_origin_parameters, radial_large_components,radial_small_components
        use radial_orbital_interpolation, only: intrpq
        use radial_orbital_orthogonalization, only: orthsc
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: wavefunction_input_unit
        integer, intent(out) :: load_status

        ! Local declarations
        integer :: orbital_index, radial_index, abs_orbital_kappa, loaded_orbital_count, read_status, file_principal_n, &
            & file_kappa, file_grid_size
        integer :: ierr
        real(real64) :: charge_times_alpha_squared, kappa_squared, file_energy_parameter, file_origin_coefficient, &
            & original_norm_squared
        real(real64), pointer :: file_large_component(:), file_small_component(:), file_radial_grid(:)
        write (output_unit, *) "Loading Radial WaveFunction File ..."
        call alloc(radial_large_components, max_radial_points, orbital_count, "PF", "LODRWF")
        call alloc(radial_small_components, max_radial_points, orbital_count, "QF", "LODRWF")

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
                    stop "lodrwf: Inappropriate gamma"
                end if
            end if
        end do
        if (ldbpr(3)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,"From SUBROUTINE LODRWF:"/," Orbital",8x,"Eigenvalue",19x,"Norm")')
        end if
        loaded_orbital_count = 0
        ! Load each binary record in file order and install it only in matching orbital slots.
        do
            read (wavefunction_input_unit, iostat=read_status) file_principal_n, file_kappa, file_energy_parameter, file_grid_size
            if (read_status /= 0) exit
            call alloc (file_large_component,file_grid_size, "PA", "LODRWF")
            call alloc (file_small_component,file_grid_size, "QA", "LODRWF")
            call alloc (file_radial_grid,file_grid_size, "RA", "LODRWF")

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
                call intrpq(file_large_component, file_small_component, file_grid_size, file_radial_grid, orbital_index, &
                    & original_norm_squared)
                if (ldbpr(3)) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)') &
                        orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), &
                        & orbital_energy_parameters(orbital_index), original_norm_squared
                end if
                loaded_orbital_count = loaded_orbital_count + 1
            end do
            call dalloc (file_large_component, "PA", "LODRWF" )
            call dalloc (file_small_component, "QA", "LODRWF")
            call dalloc (file_radial_grid, "RA", "LODRWF")
        end do
        if (ldbpr(3)) then
            call ensure_debug_output()
            write (angular_debug_unit, *) " orbitals renormalised;"
        end if
        if (loaded_orbital_count < orbital_count) then
            load_status = 1
            return
        end if
        call orthsc
        if (ldbpr(3)) then
            call ensure_debug_output()
            write (angular_debug_unit, *) " orbitals orthogonalised and renormalised;"
        end if
        load_status = 0
    end subroutine lodrwf
end module radial_wavefunction_records
