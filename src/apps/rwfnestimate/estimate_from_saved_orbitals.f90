! Read existing radial orbitals as initial estimates on the current grid.
module estimate_from_saved_orbitals_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: estimate_from_saved_orbitals
contains

    ! This subroutine loads  radial wavefunctions from the  .rwf  file
    ! and performs some related setup.
    !
    ! Written by Farid A. Parpia. Last revision: 18 Dec 1992
    subroutine estimate_from_saved_orbitals(matched_orbital_indices, matched_subshell_count, file_path)
        ! Module dependencies
        use radial_wavefunction_format, only: wavefunction_layout, read_wavefunction_header, read_orbital_metadata, &
            open_wavefunction_input, read_orbital_payload
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: alloc, dalloc
        use debug_output_settings, only: rmcdhf_debug_flags
        use initial_orbital_selection, only: max_orbitals, orbital_initialized
        use orbital_csf_basis, only: orbital_energy_parameters, orbital_principal_n, orbital_kappa, orbital_symmetry_labels
        use radial_wavefunctions, only: orbital_origin_parameters
        use initial_orbital_source, only: orbital_source_description
        use radial_orbital_interpolation, only: interpolate_radial_orbital
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: wavefunction_input_unit
        integer , intent(in) :: matched_subshell_count
        character, intent(in) :: file_path*(*)
        integer , intent(in) :: matched_orbital_indices(max_orbitals)
        integer :: open_status, io_status, saved_principal_n, saved_kappa, saved_grid_point_count, matched_subshell_index, &
            & orbital_index, file_path_length
        real(real64) :: saved_orbital_energy, radial_norm
        real(real64), pointer :: saved_large_component(:), saved_small_component(:), saved_radial_grid(:)
        type(wavefunction_layout) :: layout

        call open_wavefunction_input(wavefunction_input_unit, file_path, open_status)
        if (open_status /= 0) error stop "Cannot open radial wavefunction: "//trim(file_path)
        ! Check the file; if not as expected, try again
        call read_wavefunction_header(wavefunction_input_unit, layout, io_status)
        if (io_status /= 0) error stop "Invalid radial wavefunction header"
        ! Read orbital information from Read Orbitals File; write summary
        ! to  .dbg  file if option set
        if (rmcdhf_debug_flags(3)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,"From SUBROUTINE FRMRWF:"/," Orbital",8x,"Eigenvalue",19x,"Norm")')
        end if
        do
            call read_orbital_metadata(layout, saved_principal_n, saved_kappa, &
                saved_orbital_energy, saved_grid_point_count, io_status)
            if (is_iostat_end(io_status)) exit
            if (io_status /= 0) error stop "Invalid orbital metadata"
            do matched_subshell_index = 1, matched_subshell_count
                orbital_index = matched_orbital_indices(matched_subshell_index)
                if (.not.(.not.orbital_initialized(orbital_index) .and. &
                    orbital_principal_n(orbital_index)==saved_principal_n .and. &
                    & orbital_kappa(orbital_index)==saved_kappa)) &
                    cycle
                orbital_energy_parameters(orbital_index) = saved_orbital_energy
                call alloc (saved_large_component, saved_grid_point_count, "saved_large_component", "estimate_from_saved_orbitals")
                call alloc (saved_small_component, saved_grid_point_count, "saved_small_component", "estimate_from_saved_orbitals")
                call alloc (saved_radial_grid, saved_grid_point_count, "saved_radial_grid", "estimate_from_saved_orbitals")
                call read_orbital_payload(layout, orbital_origin_parameters(orbital_index), &
                    saved_large_component, saved_small_component, saved_radial_grid)
                call interpolate_radial_orbital (saved_large_component, saved_small_component, saved_grid_point_count, &
                    & saved_radial_grid, orbital_index, radial_norm)
                if (rmcdhf_debug_flags(3)) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)') &
                        orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), &
                        & orbital_energy_parameters(orbital_index), radial_norm
                end if
                call dalloc (saved_large_component, "saved_large_component", "estimate_from_saved_orbitals")
                call dalloc (saved_small_component, "saved_small_component", "estimate_from_saved_orbitals")
                call dalloc (saved_radial_grid, "saved_radial_grid", "estimate_from_saved_orbitals")
                file_path_length = len_trim(file_path)
                orbital_initialized(orbital_index) = .true.
                orbital_source_description(orbital_index) = file_path(1:3)
                exit
            end do
        end do
        if (rmcdhf_debug_flags(3)) then
            call ensure_debug_output()
            write (angular_debug_unit, *) " orbitals renormalised;"
        end if
        close(wavefunction_input_unit)
        return
        return
    end subroutine estimate_from_saved_orbitals
end module estimate_from_saved_orbitals_interface
