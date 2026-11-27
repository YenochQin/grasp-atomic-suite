! Transform the initial and final wavefunctions to biorthogonal orbital sets for transition
! calculations.
! Counter-transform the mixing coefficients so the initial and final atomic states are preserved.
! General references:
!
! P.A. Malmqvist, Int.J. of Quantum Chemistry, XXX, 479-94 (1986)
! J. Olsen, M.R. Godefroid, P. Jonsson, P.A. Malmqvist and
! C. Froese Fischer, Phys. Rev. E, 4499 (1995)
! Program written by
!
! Per Jonsson, Department of Computer Science
!              Vanderbilt University, USA
!
! Date June 1996
!
! Modified by Gediminas Gaigalas for new spin-angular integration
! and for reducing usage of CPU memory.        NIST, October 2017
program biorthtransform_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: int64, output_unit, real64
    use grasp_parameters,   only: max_orbitals, max_radial_points
    use calculation_input_defaults,       only: nondefault_options_flag, angular_output_dump_flag
    use biorthogonal_subshell_layout,         only: max_biorthogonal_kappa_groups, kappa_group_count, &
        initial_orbital_counts_by_group, final_orbital_counts_by_group
    use biorthogonal_orbital_sets,         only: initial_orbital_count, final_orbital_count
    use radial_wavefunctions,          only: initial_radial_large_components, initial_radial_small_components, &
        final_radial_large_components, final_radial_small_components
    use biorthogonal_orbital_rotations,         only: initial_orbital_rotation_values, final_orbital_rotation_values
    use biorthogonal_state_files, only: biorthogonal_pair, files_for_biorthogonal_pair, &
        begin_biorthogonal_output, complete_biorthogonal_output
    use csf_descriptor_reader, only: loaded_csf_basis_id
    use yes_no_input, only: read_yes_no
    use machine_constants, only: initialize_machine_constants
    use physical_constants, only: initialize_physical_constants
    use isotope_data_input, only: read_isotope_data

    ! Procedure interfaces
    use configure_radial_grid_interface, only: configure_radial_grid

    ! Module dependencies
    use radial_grid_setup, only: initialize_radial_grid
    use radial_stencil_weights, only: initialize_radial_stencil_weights

    ! Procedure interfaces
    use calculation_state_files, only: read_state_pair
    use select_biorthogonal_symmetries_interface, only: select_biorthogonal_symmetries
    use load_merged_csf_basis_interface, only: load_merged_csf_basis
    use copy_state_csf_data_interface, only: copy_state_csf_data
    use prepare_biorthogonal_subshells_interface, only: prepare_biorthogonal_subshells
    use read_initial_radial_wavefunctions_interface, only: read_initial_radial_wavefunctions
    use read_final_radial_wavefunctions_interface, only: read_final_radial_wavefunctions
    use calculate_orbital_overlap_matrix_interface, only: calculate_orbital_overlap_matrix
    use calculate_subshell_overlap_interface, only: calculate_subshell_overlap
    use transform_orbital_sets_interface, only: transform_orbital_sets
    use write_transformed_radial_wavefunctions_interface, only: write_transformed_radial_wavefunctions
    use prepare_biorthogonal_mcp_interface, only: prepare_biorthogonal_mcp

    ! Module dependencies
    use run_timer, only: begin_run_timer, finish_run_timer
    implicit none (type, external)
    integer, parameter :: orbital_transformation_workspace_size = 10000000
    integer :: ninshli(max_biorthogonal_kappa_groups), ninshlf(max_biorthogonal_kappa_groups)
    integer :: global_trace_level, local_trace_level, trace_level, mixing_file_kind, initial_core_subshell_count, &
        & final_core_subshell_count, mxl
    integer(int64) :: start_clock_count
    real(real64) :: orbital_transformation_workspace(orbital_transformation_workspace_size)
    real(real64) :: cishl(max_biorthogonal_kappa_groups*max_biorthogonal_kappa_groups)
    real(real64)   :: s(max_orbitals, max_orbitals)
    real(real64) :: cfshl(max_biorthogonal_kappa_groups*max_biorthogonal_kappa_groups)
    logical :: yes
    character :: name(2)*24
    character(len=128) :: isotope_file_path

    type(biorthogonal_pair) :: pair
    character(64) :: basis_ids(2)

    call begin_run_timer (start_clock_count, "BIORTHTRANSFORM")
    isotope_file_path = "isodata"
    ! Debug flags
    global_trace_level = 0
    local_trace_level = 0
    trace_level = max0(local_trace_level,global_trace_level)
    write(*,*) "This program transforms the initial and final wave"
    write(*,*) "functions so that standard tensor albegra can be"
    write(*,*) "used in evaluation of the transition parameters"
    write(*,*) "Input files:  isodata, name1_desc.parquet, name1.wfn, name1.mix/name1.cimix"
    write(*,*) "              name2_desc.parquet, name2.wfn, name2.mix/name2.cimix         "
    print*,"              name1.biorth.angular.cache, name2.biorth.angular.cache (optional angular files)"
    write(*,*) "Output: biorth/<initial>/<final>/<mcdhf|ci>/ (two orbitals, two mixing files, manifest)"
    write(*,*) "              name1.biorth.angular.cache, name2.biorth.angular.cache (angular files)"
    write(*,*)

    write (output_unit, *) "Default settings?"
    yes = read_yes_no()
    write (output_unit, *)
    if (yes) then
        nondefault_options_flag = 0
        angular_output_dump_flag = 1
    else
        nondefault_options_flag = 1
        write (output_unit, *) "Dump angular data on file?"
        yes = read_yes_no()
        write (output_unit, *)
        if (yes) then
            angular_output_dump_flag = 1
        else
            angular_output_dump_flag = 0
        end if
    end if
    write (output_unit, *) "Input from a CI calculation?"
    yes = read_yes_no()
    write (output_unit, *)
    if (yes) then
        mixing_file_kind = 0
    else
        mixing_file_kind = 1
    end if
    ! Initialize machine-dependent constants.
    call initialize_machine_constants
    ! Initialize physical constants.
    call initialize_physical_constants
    ! Read the nuclear parameters from isodata.
    call read_isotope_data (isotope_file_path)
    ! Choose the radial grid parameters.
    call configure_radial_grid
    ! Build the radial grid.
    call initialize_radial_grid
    ! Initialize the radial integration and derivative stencils.
    call initialize_radial_stencil_weights
    ! Resolve the input paths and prepare output for the transformed orbitals and mixing coefficients.
    call read_state_pair(name(1), name(2))
    if (name(1) == name(2)) then
        write(*, "(a)") "Same state on both sides. Run rtransition with the original orbital and mixing files."
        call finish_run_timer (start_clock_count, "BIORTHTRANSFORM")
        stop
    end if
    pair = files_for_biorthogonal_pair(name, mixing_file_kind == 0, ".")
    call select_biorthogonal_symmetries()
    ! Load the initial-state CSF descriptors.
    call load_merged_csf_basis (name(1), initial_core_subshell_count,1)
    basis_ids(1) = loaded_csf_basis_id
    ! Preserve the initial-state data before loading the final state.
    call copy_state_csf_data (1)
    ! Load the final-state CSF descriptors.
    call load_merged_csf_basis (name(2), final_core_subshell_count,2)
    basis_ids(2) = loaded_csf_basis_id
    ! Preserve the final-state data for the transformation.
    call copy_state_csf_data (2)
    ! Determine the number of kappa quantum numbers and
    ! the number of orbitals for each kappa quantum number
    ! for the initial state and final states
    call prepare_biorthogonal_subshells (global_trace_level, initial_core_subshell_count, final_core_subshell_count)
    call begin_biorthogonal_output(pair)
    ! Read the the radial orbitals for the initial state
    call read_initial_radial_wavefunctions (name(1), global_trace_level)
    ! Read the the radial orbitals for the initial state
    call read_final_radial_wavefunctions (name(2), global_trace_level)
    ! Calculate the radial overlap matrices
    write (*, *)
    write (*, *) " ******************************************"
    write (*, *) "  Overlap matrix before orbital rotations"
    write (*, *) " *****************************************"
    write (*, *)

    call calculate_orbital_overlap_matrix

    call calculate_subshell_overlap (s, initial_orbital_count, final_orbital_count)

    ! Once we have the overlap matrices
    ! we can manipulate the initial and final state separately.
    mxl = kappa_group_count
    ! . Calculate biorthonormal orbitals, and orbital matrix
    ! . for counter transformation of CI coefficients.
    call transform_orbital_sets (initial_radial_large_components, initial_radial_small_components, &
        initial_orbital_counts_by_group, ninshli, &
        final_radial_large_components, final_radial_small_components, &
        final_orbital_counts_by_group, ninshlf, max_radial_points, kappa_group_count, orbital_transformation_workspace, &
            orbital_transformation_workspace_size, &
            & global_trace_level, &
        cishl, initial_orbital_rotation_values, cfshl, final_orbital_rotation_values)
    write (*, *)
    write (*, *) " ****************************************"
    write (*, *) "  Overlap matrix after orbital rotations"
    write (*, *) " ****************************************"
    write (*, *)

    call calculate_orbital_overlap_matrix
    ! Write the transformed radial functions to file
    call write_transformed_radial_wavefunctions (pair%sides)
    ! Obtain one-electron coupling coefficients for the initial state.
    ! The coefficients are dumped on files one kappa in turn and
    ! thus the different kappa can be manipulated independently.
    ! The interface with the transformation part is in the routine mcp
    call prepare_biorthogonal_mcp (name(1), 1, global_trace_level, mixing_file_kind, pair%sides(1)%mixing)
    ! Obtain one-electron coupling coefficients for the final state.
    ! The coefficients are dumped on files one kappa in turn and
    ! thus the different kappa can be manipulated independently.
    ! The interface with the transformation part is in the routine mcp
    call prepare_biorthogonal_mcp (name(2), 2, global_trace_level, mixing_file_kind, pair%sides(2)%mixing)

    call complete_biorthogonal_output(pair, name, basis_ids)
    call finish_run_timer (start_clock_count, "BIORTHTRANSFORM")

    stop
end program biorthtransform_main
