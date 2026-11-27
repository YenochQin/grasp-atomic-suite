! Handle records shared by the serial and distributed RCI restart files.
module rci_restart_records
    use, intrinsic :: iso_fortran_env, only: iostat_end, real64
    use input_error_handling, only: abort_input, fail_input_io
    use csf_descriptor_reader, only: loaded_csf_basis_id
    implicit none (type, external)

    private
    public :: open_rci_restart_file, write_rci_restart_identity, read_rci_restart_header, scan_rci_matrix_block

    character(*), parameter :: restart_signature = "R92RS2"
    character(*), parameter :: digest_digits = "0123456789abcdef"

contains

    subroutine open_rci_restart_file(unit, path, is_restart, abort_handler)
        integer, intent(out) :: unit
        character(*), intent(in) :: path
        logical, intent(in) :: is_restart
        procedure(abort_input), optional :: abort_handler

        integer :: status
        character(len(restart_signature)) :: signature

        if (is_restart) then
            open(newunit=unit, file=path, status="old", form="unformatted", position="rewind", iostat=status)
            if (status /= 0) call fail_input_io(status, -1, "RCI cannot open restart file "//path, abort_handler)
            read(unit, iostat=status) signature
            if (status /= 0) call fail_input_io(status, unit, "RCI restart signature", abort_handler)
            if (signature == "R92RES") call fail_input_io(0, unit, &
                "Legacy R92RES restart has no CSF basis digest. Start a new RCI calculation.", abort_handler)
            if (signature /= restart_signature) call fail_input_io(0, unit, &
                "RCI restart signature must be "//restart_signature, abort_handler)
        else
            open(newunit=unit, file=path, status="unknown", form="unformatted", position="rewind", iostat=status)
            if (status /= 0) call fail_input_io(status, -1, "RCI cannot open matrix file "//path, abort_handler)
        end if
    end subroutine open_rci_restart_file

    subroutine write_rci_restart_identity(unit, abort_handler)
        use atomic_calculation_state, only: electron_count
        use orbital_csf_basis, only: csf_count, orbital_count
        use hamiltonian_symmetry_blocks, only: symmetry_block_count

        integer, intent(in) :: unit
        procedure(abort_input), optional :: abort_handler
        integer :: status

        if (verify(loaded_csf_basis_id, digest_digits) /= 0) &
            call fail_input_io(0, unit, "RCI cannot write restart without a valid CSF basis digest", abort_handler)
        write(unit, iostat=status) restart_signature
        if (status /= 0) call fail_input_io(status, unit, "RCI writing restart signature", abort_handler)
        write(unit, iostat=status) electron_count, csf_count, orbital_count, symmetry_block_count, loaded_csf_basis_id
        if (status /= 0) call fail_input_io(status, unit, "RCI writing restart basis identity", abort_handler)
    end subroutine write_rci_restart_identity

    ! The original header and radial records were written by Farid A. Parpia and Xinghong He.
    subroutine read_rci_restart_header(unit, abort_handler)
        use grasp_parameters, only: max_radial_points, radial_workspace_capacity
        use hamiltonian_correction_flags, only: use_first_order_csf_perturbation, include_normal_mass_shift, &
            include_specific_mass_shift, include_transverse_interaction, include_vacuum_polarization
        use atomic_calculation_state, only: left_electron_count, electron_count, nuclear_charge, &
            nuclear_mass_in_electron_masses, speed_of_light_au, working_tolerance
        use radial_grid_state, only: radial_point_count, radial_grid, radial_grid_derivative, radial_derivative_over_radius, &
            radial_grid_scale, transformed_grid_step, asymptotic_grid_step
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
        use nuclear_potential_state, only: negative_radius_times_potential, first_grid_index_beyond_fermi_c
        use orbital_csf_basis, only: csf_count, orbital_count, orbital_energy_parameters, orbital_origin_exponents
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters, radial_large_components, &
            radial_small_components
        use breit_frequency_scale, only: transverse_photon_frequency_scale
        use hamiltonian_symmetry_blocks, only: symmetry_block_count, csf_counts_by_block
        use csf_block_cutoffs, only: csf_cutoff_by_block
        use radial_stencil_weights, only: initialize_radial_stencil_weights

        integer, intent(in) :: unit
        procedure(abort_input), optional :: abort_handler

        integer :: status, file_csfs, file_orbitals, file_blocks, record_count, padded_points, orbital, last_point
        character(len(loaded_csf_basis_id)) :: file_basis_id
        character(16) :: orbital_label

        read(unit, iostat=status) left_electron_count, file_csfs, file_orbitals, file_blocks, file_basis_id
        if (status /= 0) call fail_input_io(status, unit, "RCI restart basis identity", abort_handler)
        if (left_electron_count /= electron_count .or. file_csfs /= csf_count .or. file_orbitals /= orbital_count .or. &
            file_blocks /= symmetry_block_count) &
            call fail_input_io(0, unit, "RCI restart basis counts differ from the current basis", abort_handler)
        if (verify(file_basis_id, digest_digits) /= 0) &
            call fail_input_io(0, unit, "RCI restart has an invalid CSF basis digest", abort_handler)
        if (file_basis_id /= loaded_csf_basis_id) &
            call fail_input_io(0, unit, "RCI restart basis digest differs from the current CSF basis", abort_handler)

        read(unit, iostat=status) nuclear_charge, nuclear_mass_in_electron_masses
        if (status /= 0) call fail_input_io(status, unit, "RCI restart nuclear charge and mass", abort_handler)
        read(unit, iostat=status) nuclear_model_parameter_count
        if (status /= 0) call fail_input_io(status, unit, "RCI restart nuclear model parameter count", abort_handler)
        if (nuclear_model_parameter_count /= 0 .and. nuclear_model_parameter_count /= size(nuclear_model_parameters)) &
            call fail_input_io(0, unit, "RCI restart invalid nuclear model parameter count", abort_handler)
        backspace(unit, iostat=status)
        if (status /= 0) call fail_input_io(status, unit, "RCI restart positioning nuclear model record", abort_handler)
        read(unit, iostat=status) record_count, nuclear_model_parameters(:nuclear_model_parameter_count)
        if (status /= 0) call fail_input_io(status, unit, "RCI restart nuclear model parameters", abort_handler)
        if (record_count /= nuclear_model_parameter_count) &
            call fail_input_io(0, unit, "RCI restart nuclear model parameter count changed on reread", abort_handler)

        read(unit, iostat=status) radial_point_count
        if (status /= 0) call fail_input_io(status, unit, "RCI restart radial point count", abort_handler)
        if (radial_point_count < 1 .or. radial_point_count > max_radial_points) &
            call fail_input_io(0, unit, "RCI restart invalid radial point count", abort_handler)
        backspace(unit, iostat=status)
        if (status /= 0) call fail_input_io(status, unit, "RCI restart positioning nuclear potential record", abort_handler)
        read(unit, iostat=status) record_count, negative_radius_times_potential(:radial_point_count), &
            first_grid_index_beyond_fermi_c
        if (status /= 0) call fail_input_io(status, unit, "RCI restart nuclear potential", abort_handler)
        if (record_count /= radial_point_count) &
            call fail_input_io(0, unit, "RCI restart radial point count changed on reread", abort_handler)

        read(unit, iostat=status) speed_of_light_au, use_first_order_csf_perturbation, csf_cutoff_by_block(:symmetry_block_count), &
            include_transverse_interaction, transverse_photon_frequency_scale, include_vacuum_polarization, &
            include_normal_mass_shift, include_specific_mass_shift
        if (status /= 0) call fail_input_io(status, unit, "RCI restart Hamiltonian options", abort_handler)
        if (any(csf_cutoff_by_block(:symmetry_block_count) < 0) .or. &
            any(csf_cutoff_by_block(:symmetry_block_count) > csf_counts_by_block(1:symmetry_block_count))) &
            call fail_input_io(0, unit, "RCI restart CSF cutoff is outside its symmetry block", abort_handler)

        padded_points = radial_point_count + radial_workspace_capacity - max_radial_points
        read(unit, iostat=status) radial_grid_scale, transformed_grid_step, asymptotic_grid_step, radial_grid(:padded_points), &
            radial_grid_derivative(:padded_points), radial_derivative_over_radius(:padded_points)
        if (status /= 0) call fail_input_io(status, unit, "RCI restart radial grid", abort_handler)
        working_tolerance = transformed_grid_step**6
        call initialize_radial_stencil_weights

        allocate(radial_large_components(max_radial_points, orbital_count), &
            radial_small_components(max_radial_points, orbital_count), stat=status)
        if (status /= 0) call fail_input_io(status, unit, "RCI restart allocating radial orbitals", abort_handler)
        do orbital = 1, orbital_count
            write(orbital_label, "(i0)") orbital
            read(unit, iostat=status) orbital_energy_parameters(orbital), orbital_origin_exponents(orbital), &
                orbital_origin_parameters(orbital), last_point
            if (status /= 0) call fail_input_io(status, unit, "RCI restart metadata for orbital "//trim(orbital_label), &
                abort_handler)
            if (last_point < 1 .or. last_point > radial_point_count) &
                call fail_input_io(0, unit, "RCI restart invalid radial point count for orbital "//trim(orbital_label), &
                    abort_handler)
            orbital_last_grid_index(orbital) = last_point
            read(unit, iostat=status) radial_large_components(:last_point, orbital), radial_small_components(:last_point, orbital)
            if (status /= 0) call fail_input_io(status, unit, "RCI restart radial components for orbital "//trim(orbital_label), &
                abort_handler)
            ! Radial differentiation uses the full active grid, including the unsaved tail.
            radial_large_components(last_point+1:radial_point_count, orbital) = 0.0_real64
            radial_small_components(last_point+1:radial_point_count, orbital) = 0.0_real64
        end do
    end subroutine read_rci_restart_header

    ! A clean end of file at a record boundary is a partial checkpoint and can be resumed.
    subroutine scan_rci_matrix_block(unit, block, csfs, cutoff, rank, ranks, header_status, rows, last_shift, abort_handler)
        use sparse_hamiltonian_storage, only: stored_matrix_element_count
        use eigenstate_expansions, only: hamiltonian_energy_shift
        integer, intent(in) :: unit, block, csfs, cutoff, rank, ranks
        integer, intent(out) :: header_status, rows
        real(real64), intent(inout) :: last_shift
        procedure(abort_input), optional :: abort_handler

        integer :: status, file_csfs, file_cutoff, file_rank, file_ranks, row, entries, record_count, column
        integer, allocatable :: columns(:)
        real(real64) :: unused_value, diagonal, shift
        character(72) :: context

        rows = 0
        read(unit, iostat=header_status) file_csfs, file_cutoff, file_rank, file_ranks
        if (header_status == iostat_end) return
        write(context, "(a,i0)") "RCI restart matrix block header for block ", block
        if (header_status /= 0) call fail_input_io(header_status, unit, context, abort_handler)
        if (file_csfs /= csfs .or. file_cutoff /= cutoff .or. file_rank /= rank .or. file_ranks /= ranks) &
            call fail_input_io(0, unit, trim(context)//" does not match CSF cutoff or rank layout", abort_handler)

        allocate(columns(csfs), stat=status)
        if (status /= 0) call fail_input_io(status, unit, "RCI restart allocating matrix column indices", abort_handler)
        do row = rank + 1, csfs, ranks
            read(unit, iostat=status) entries
            if (status == iostat_end) exit
            write(context, "(a,i0,a,i0)") "RCI restart matrix row entry count, block ", block, ", row ", row
            if (status /= 0) call fail_input_io(status, unit, context, abort_handler)
            if (entries < 1 .or. entries > row) call fail_input_io(0, unit, context, abort_handler)
            backspace(unit, iostat=status)
            if (status /= 0) call fail_input_io(status, unit, "RCI restart positioning "//trim(context), abort_handler)
            read(unit, iostat=status) record_count, shift, (unused_value, column=2,entries), diagonal, columns(:entries)
            write(context, "(a,i0,a,i0)") "RCI restart matrix row, block ", block, ", row ", row
            if (status /= 0) call fail_input_io(status, unit, context, abort_handler)
            if (record_count /= entries) &
                call fail_input_io(0, unit, trim(context)//" entry count changed on reread", abort_handler)
            if (columns(1) < 1 .or. columns(entries) /= row .or. any(columns(2:entries) <= columns(:entries-1))) &
                call fail_input_io(0, unit, trim(context)//" has invalid column indices", abort_handler)
            rows = rows + 1
            last_shift = shift
            stored_matrix_element_count = stored_matrix_element_count + entries
            hamiltonian_energy_shift = hamiltonian_energy_shift + diagonal
        end do
    end subroutine scan_rci_matrix_block
end module rci_restart_records
