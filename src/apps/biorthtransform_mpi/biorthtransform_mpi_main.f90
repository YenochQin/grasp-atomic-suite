! Coordinate the MPI calculation of biorthogonal orbitals and transformed mixing coefficients.
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
program biorthtransform_mpi_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use grasp_parameters,   only: max_orbitals, max_radial_points
    use calculation_input_defaults,       only: nondefault_options_flag, angular_output_dump_flag
    use biorthogonal_subshell_layout,         only: max_biorthogonal_kappa_groups, kappa_group_count, &
        initial_orbital_counts_by_group, final_orbital_counts_by_group
    use biorthogonal_orbital_sets,         only: initial_orbital_count, final_orbital_count
    use radial_wavefunctions,          only: initial_radial_large_components, initial_radial_small_components, &
        final_radial_large_components, final_radial_small_components
    use biorthogonal_orbital_rotations,         only: initial_orbital_rotation_values, final_orbital_rotation_values
    use mpi_state, only: mpi_barrier, mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_process_count, mpi_rank, &
        processor_name, processor_name_length
    use biorthogonal_state_files, only: biorthogonal_pair, files_for_biorthogonal_pair, &
        begin_biorthogonal_output, complete_biorthogonal_output
    use csf_descriptor_reader, only: loaded_csf_basis_id
    use mpi_io_errors, only: abort_input_mpi
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
    use mpi_calculation_state, only: read_mpi_state_pair
    use select_biorthogonal_symmetries_interface, only: select_biorthogonal_symmetries
    use load_merged_csf_basis_interface, only: load_merged_csf_basis
    use copy_state_csf_data_interface, only: copy_state_csf_data
    use biorthogonal_transition_selection, only: ntrans, itjqj1, itjqj2, jqj1, jqj2
    use prepare_biorthogonal_subshells_interface, only: prepare_biorthogonal_subshells
    use read_initial_radial_wavefunctions_interface, only: read_initial_radial_wavefunctions
    use read_final_radial_wavefunctions_interface, only: read_final_radial_wavefunctions
    use calculate_orbital_overlap_matrix_interface, only: calculate_orbital_overlap_matrix
    use calculate_subshell_overlap_interface, only: calculate_subshell_overlap
    use transform_orbital_sets_interface, only: transform_orbital_sets
    use write_transformed_radial_wavefunctions_interface, only: write_transformed_radial_wavefunctions

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table

    ! Procedure interfaces
    use read_biorthogonal_mcp_interface, only: read_biorthogonal_mcp

    ! Module dependencies
    use mpi_utilities, only: initialize_mpi_runtime_with_status
    use mpi_utilities, only: finish_mpi_application

    ! Procedure interfaces
    use load_biorthogonal_csf_basis_mpi_interface, only: load_biorthogonal_csf_basis_mpi

    ! Module dependencies
    use mpi_change_directory, only: change_mpi_directory

    ! Procedure interfaces
    use write_biorthogonal_mcp_interface, only: write_biorthogonal_mcp
    implicit none (type, external)
    integer, parameter :: orbital_transformation_workspace_size = 10000000
    integer :: ninshli(max_biorthogonal_kappa_groups), ninshlf(max_biorthogonal_kappa_groups)
    integer :: ntestg, ntestl, ntest, mixing_file_kind, initial_core_subshell_count, final_core_subshell_count, mxl
    integer :: permanent_directory_length,scratch_directory_length,iii,unused_core_subshell_count
    integer(int64) :: start_clock_count
    real(real64) :: orbital_transformation_workspace(orbital_transformation_workspace_size)
    real(real64) :: cishl(max_biorthogonal_kappa_groups*max_biorthogonal_kappa_groups)
    real(real64)   :: s(max_orbitals, max_orbitals)
    real(real64) :: cfshl(max_biorthogonal_kappa_groups*max_biorthogonal_kappa_groups)
    logical :: yes
    character :: name(2)*24
    character(len=3)  ::  block_label_text
    character(len=128) :: isotope_file_path
    character(len=128) :: fullname(2)
    character (len = 128) :: startup_directory,permdir,tmpdir
    ! Initialize MPI and set up the I/O and scratch directories.
    type(biorthogonal_pair) :: pair
    character(64) :: basis_ids(2)

    call initialize_mpi_runtime_with_status (mpi_rank, mpi_process_count, processor_name, processor_name_length, &
        & start_clock_count,            &
        startup_directory, permdir, tmpdir, "BIORTHTRANSFORM_MPI")
    write (block_label_text, "(I3.3)") mpi_rank
    print*, tmpdir , " = tmpdir"

    isotope_file_path = trim(startup_directory)//"/isodata"
    permanent_directory_length = len_trim (permdir)
    scratch_directory_length = len_trim (tmpdir)

    ! Debug flags
    ntestg = 0
    ntestl = 0
    ntest = max0(ntestl,ntestg)
    if (mpi_rank == 0) then
        print *
        print *

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
    end if !myid=0
    call mpi_bcast (nondefault_options_flag,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
    call mpi_bcast (angular_output_dump_flag,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
    call mpi_bcast (mixing_file_kind,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
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
    call read_mpi_state_pair(name)
    if (name(1) == name(2)) then
        if (mpi_rank == 0) write(*, "(a)") "Same state on both sides. Run rtransition with the original orbital and mixing files."
        call finish_mpi_application (mpi_rank, start_clock_count, "BIORTHTRANSFORM_MPI")
        stop
    end if
    pair = files_for_biorthogonal_pair(name, mixing_file_kind == 0, trim(startup_directory))
    if (mpi_rank == 0) call select_biorthogonal_symmetries()
    call mpi_bcast(ntrans, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
    if (ntrans == 1) then
        call mpi_bcast(jqj1, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
        call mpi_bcast(itjqj1, jqj1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
        call mpi_bcast(jqj2, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
        call mpi_bcast(itjqj2, jqj2, mpi_integer, 0, mpi_comm_world, mpi_error_code)
    end if
    print *, "DBG: myid=",mpi_rank, " name(1)=", name(1)
    print *, "DBG: myid=",mpi_rank, " name(2)=", name(2)
    fullname(1) = trim(startup_directory) // "/" // trim(name(1))
    fullname(2) = trim(startup_directory) // "/" // trim(name(2))
    print *, "DBG: myid=",mpi_rank, " startdir =", startup_directory
    print *, "DBG: myid=",mpi_rank, " fullname(1)=", fullname(1)
    print *, "DBG: myid=",mpi_rank, " fullname(2)=", fullname(2)
    ! Load the initial-state CSF descriptors.
    call load_merged_csf_basis (fullname(1), initial_core_subshell_count,1)
    basis_ids(1) = loaded_csf_basis_id
    ! Preserve the initial-state data before loading the final state.
    call copy_state_csf_data (1)
    ! Load the final-state CSF descriptors.
    call load_merged_csf_basis (fullname(2), final_core_subshell_count,2)
    basis_ids(2) = loaded_csf_basis_id
    ! Preserve the final-state data for the transformation.
    call copy_state_csf_data (2)
    ! Determine the number of kappa quantum numbers and
    ! the number of orbitals for each kappa quantum number
    ! for the initial state and final states
    call prepare_biorthogonal_subshells (ntestg, initial_core_subshell_count, final_core_subshell_count)
    if (mpi_rank == 0) call begin_biorthogonal_output(pair, abort_input_mpi)
    call mpi_barrier(mpi_comm_world, mpi_error_code)
    ! Read the the radial orbitals for the initial state
    call read_initial_radial_wavefunctions (fullname(1), ntestg)
    ! Read the the radial orbitals for the initial state
    call read_final_radial_wavefunctions (fullname(2), ntestg)
    ! Calculate the radial overlap matrices
    if (mpi_rank == 0) then
        write (*, *)
        write (*, *) " ******************************************"
        write (*, *) "  Overlap matrix before orbital rotations"
        write (*, *) " *****************************************"
        write (*, *)
    end if!myid=0

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
            & ntestg, &
        cishl, initial_orbital_rotation_values, cfshl, final_orbital_rotation_values)
    if (mpi_rank == 0) then
        write (*, *)
        write (*, *) " ****************************************"
        write (*, *) "  Overlap matrix after orbital rotations"
        write (*, *) " ****************************************"
        write (*, *)
    end if !myid=0

    call calculate_orbital_overlap_matrix
    ! Write the transformed radial functions to file
    if (mpi_rank == 0) call write_transformed_radial_wavefunctions(pair%sides)
    ! Obtain one-electron coupling coefficients for the initial state.
    ! The coefficients are dumped on files one kappa in turn and
    ! thus the different kappa can be manipulated independently.
    ! The interface with the transformation part is in the routine mcp
    iii = len_trim(startup_directory)
    call change_mpi_directory(trim(startup_directory),iii,mpi_error_code)
    call load_biorthogonal_csf_basis_mpi(name(1), unused_core_subshell_count,3)
    iii = len_trim(tmpdir)
    call change_mpi_directory(trim(tmpdir)//"/"//block_label_text,iii+4,mpi_error_code)
    call initialize_log_factorial_table
    call write_biorthogonal_mcp(name(1),startup_directory,1,ntestg,mixing_file_kind)
    call read_biorthogonal_mcp(name(1),startup_directory,1,ntestg,mixing_file_kind,pair%sides(1)%mixing)

    ! added by Yu Zou, Feb.18,2000 ***********
    ! R      DO K = 1,KAMAX
    ! R        close(UNIT=80+K,STATUS="DELETE")
    ! R      ENDDO
    ! added by Yu Zou, Feb.18,2000 ***********
    !
    !    Obtain one-electron coupling coefficients for the final state.
    !    The coefficients are dumped on files one kappa in turn and
    !    thus the different kappa can be manipulated independently.
    !    The interface with the transformation part is in the routine mcp
    iii = len_trim(startup_directory)
    call change_mpi_directory(trim(startup_directory),iii,mpi_error_code)
    call load_biorthogonal_csf_basis_mpi(name(2), unused_core_subshell_count,3)
    iii = len_trim(tmpdir)
    call change_mpi_directory(trim(tmpdir)//"/"//block_label_text,iii+4,mpi_error_code)
    call write_biorthogonal_mcp(name(2),startup_directory,2,ntestg,mixing_file_kind)
    call read_biorthogonal_mcp(name(2),startup_directory,2,ntestg,mixing_file_kind,pair%sides(2)%mixing)

    ! R      DO K = 1,KAMAX
    ! R        close(UNIT=80+K,STATUS="DELETE")
    ! R      ENDDO

    ! Report elapsed time.

    if (mpi_rank == 0) call complete_biorthogonal_output(pair, name, basis_ids, abort_input_mpi)
    call finish_mpi_application (mpi_rank, start_clock_count, "BIORTHTRANSFORM_MPI")
    stop
end program biorthtransform_mpi_main
