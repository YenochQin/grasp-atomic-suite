! Read the radial, orbital, and level options for the SCF calculation.
module read_scf_options_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_scf_options
contains

    !   Interactively determines the data governing the SCF problem.
    !
    ! Written by Farid A. Parpia. Last revision: 27 Dec 1992
    !
    ! Xinghong He 98-08-06
    subroutine read_scf_options(extended_optimal_level, block_labels, isotope_file_path, radial_file_path, &
        & nuclear_polarization_file_path)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: input_unit
        use calculation_file_units, only: calculation_log_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: max_orbitals
        use radial_node_count_settings, only: relative_extremum_cutoff
        use scf_damping_factors, only: orbital_damping_factors, mixing_damping_factors
        use atomic_calculation_state, only: working_tolerance, max_radial_points, selected_level_count, &
            maximum_scf_iterations, orbital_updates_per_scf_iteration, minimum_radial_solver_iterations, nuclear_charge, &
            speed_of_light_au, vacuum_speed_of_light_au
        use calculation_input_defaults, only: nondefault_options_flag
        use fixed_orbital_flags, only: is_orbital_fixed
        use standard_io_units, only: stderr_unit
        use radial_grid_state, only: transformed_grid_step, asymptotic_grid_step, radial_point_count, radial_grid_scale, radial_grid
        use orbital_inversion_flags, only: suppress_radial_inversion
        use rmcdhf_mcp_options, only: mcp_diagonal_only
        use serial_mpi_context, only: mpi_rank
        use expected_orbital_nodes, only: expected_radial_node_counts
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count
        use nuclear_potential_state, only: negative_radius_times_potential
        use nuclear_polarization_state, only: npol_message_len
        use nuclear_polarization_initialization, only: npol_init
        use orbital_csf_basis , only: orbital_kappa, orbital_symmetry_labels, orbital_principal_n, orbital_count, &
            previous_orbital_energy_changes
        use scf_orthogonality_control, only: enforce_orbital_orthogonality
        use scf_iteration_state, only: orbital_self_consistency_errors, radial_solution_methods
        use radial_wavefunctions, only: max_orbitals
        use yes_no_input, only: read_yes_no
        use isotope_data_input, only: read_isotope_data
        use radial_stencil_weights, only: initialize_radial_stencil_weights
        use radial_grid_setup, only: initialize_radial_grid
        use nuclear_coulomb_potential, only: calculate_nuclear_coulomb_potential
        use radial_wavefunction_input, only: read_radial_wavefunction_file

        ! Procedure interfaces
        use read_eol_average_options_interface, only: read_eol_average_options

        ! Module dependencies
        use integer_format, only: format_integer
        use subshell_label_selection, only: select_subshell_labels
        implicit none (type, external)

        ! Arguments
        logical, intent(out) :: extended_optimal_level
        character  :: isotope_file_path*(*)
        character  :: radial_file_path*(*)
        character(*), intent(in) :: nuclear_polarization_file_path
        character  :: block_labels(*)*8

        ! Local constants
        character(len=*), parameter :: routine_name = "read_scf_options"

        ! Local declarations
        integer :: np_status
        character(npol_message_len) :: np_message
        integer  :: indx(max_orbitals)
        integer :: i, j, iend, ibeg, text_length, nsubs, loc
        real(real64) :: odampu, cdampu
        logical :: yes
        character :: record*80, number_text*20
        ! Read the nuclear parameters from isodata.
        call read_isotope_data (isotope_file_path)
        ! Set default speed of light and grid parameters
        speed_of_light_au = vacuum_speed_of_light_au
        if (nuclear_model_parameter_count == 0) then
            radial_grid_scale = exp((-65.0_real64/16.0_real64))/nuclear_charge
            transformed_grid_step = 0.5_real64**4
            radial_point_count = min(220,max_radial_points)
        else
            !          default comes here
            ! The screening estimate could be refined to depend on nuclear charge Z.
            radial_grid_scale = 2.0e-6_real64/nuclear_charge
            transformed_grid_step = 5.0e-2_real64
            radial_point_count = max_radial_points
        end if
        asymptotic_grid_step = 0.0_real64
        ! ACCY is an estimate of the accuracy of the numerical procedures
        working_tolerance = transformed_grid_step**6
        if (nondefault_options_flag /= 0) then

            write (stderr_unit,"(A)",advance="NO")"Change the default speed of",&
                " light or radial grid parameters?  (y/n) "
            yes = read_yes_no()
            if (yes) then
                write (stderr_unit,*) "Speed of light = ",vacuum_speed_of_light_au,";", " revise ?"
                yes = read_yes_no()
                if (yes) then
                    write (stderr_unit, *) "Enter the revised value:"
                    read (input_unit, *) speed_of_light_au
                end if
                ! Choose the radial grid parameters.
                write (stderr_unit, *) "The default radial grid parameters for ",&
                    "this case are:"
                write (stderr_unit, *) " RNT = ", radial_grid_scale
                write (stderr_unit, *) " H = ", transformed_grid_step
                write (stderr_unit, *) " HP = ", asymptotic_grid_step
                write (stderr_unit, *) " N = ", radial_point_count
                write (stderr_unit, *) " revise these values?"
                yes = read_yes_no()
                if (yes) then
                    write (stderr_unit, *) "Enter RNT:"
                    read (input_unit, *) radial_grid_scale
                    write (stderr_unit, *) "Enter H:"
                    read (input_unit, *) transformed_grid_step
                    write (stderr_unit, *) "Enter HP:"
                    read (input_unit, *) asymptotic_grid_step
                    write (stderr_unit, *) "Enter N:"
                    read (input_unit, *) radial_point_count
                    ! b Revised grid
                    write (stderr_unit,*) "Revised RNT = ", radial_grid_scale
                    write (stderr_unit,*) "Revised H   = ", transformed_grid_step
                    write (stderr_unit,*) "Revised HP  = ", asymptotic_grid_step
                    write (stderr_unit,*) "Revised N   = ", radial_point_count
                end if
            end if
            ! b
            ! b read ACCY on input
            ! b
            write (stderr_unit,*) "Revise the default ACCY = ", working_tolerance
            yes = read_yes_no ()
            if (yes) then
                write (stderr_unit,*) "Enter ACCY:"
                read *, working_tolerance
                ! b Revised ACCY
                write (stderr_unit,*) "Revised ACCY = ", working_tolerance
            end if
        end if
        ! Initialize the radial integration and derivative stencils.
        call initialize_radial_stencil_weights
        ! Build the radial grid and its derivative arrays.
        call initialize_radial_grid
        ! Generate $- r \times V_ (r)$
        call calculate_nuclear_coulomb_potential
        call npol_init(nuclear_polarization_file_path, nuclear_charge, radial_point_count, radial_grid, &
            & negative_radius_times_potential, &
            & np_status, np_message)
        if (np_status /= 0) then
            write(stderr_unit,"(A)") trim(np_message)
            error stop 1
        end if
        ! Load the subshell radial wavefunction estimates
        call read_radial_wavefunction_file (radial_file_path)
        ! Set some defaults
        relative_extremum_cutoff = 0.05_real64

        ! IORDER(I) = I    ! Completely determined in GETOLD
        radial_solution_methods(:orbital_count) = 3
        suppress_radial_inversion(:orbital_count) = .true.
        ! ODAMP(:NW) = 1.D0
        orbital_damping_factors(:orbital_count) = 0.0_real64
        previous_orbital_energy_changes(:orbital_count)   = 0.0_real64
        orbital_self_consistency_errors(:orbital_count) = 0.0_real64

        where (orbital_kappa(:orbital_count) < 0)
            expected_radial_node_counts(:orbital_count) = orbital_principal_n(:orbital_count) + orbital_kappa(:orbital_count)
        else where
            expected_radial_node_counts(:orbital_count) = orbital_principal_n(:orbital_count) - orbital_kappa(:orbital_count) - 1
        end where

        if (mcp_diagonal_only) error stop &
            "RMCDHF requires ASF roots in every J/parity block. Regenerate diagonal-only MCP with rangular."
        extended_optimal_level = .true.
        call read_eol_average_options (block_labels)

        write (stderr_unit, *) "Enter the maximum number of read_scf_options cycles:"
        read (*, *) maximum_scf_iterations

        if (nondefault_options_flag==0) then
            write(calculation_log_unit,*) maximum_scf_iterations,"! Number of read_scf_options cycles"
        end if
        ! Allow the user to modify other defaults
        if (nondefault_options_flag /= 0) then
            write (stderr_unit, "(A)", advance="NO") "Modify other defaults?  (y/n) "
            yes = read_yes_no()
        else
            yes = .false.
        end if

        if (.not.yes) return
        ! From here to end, "other defaults" are handled. For simplicity
        ! We'll let node-0 do the job and then broadcast results to all
        ! nodes.
        if (mpi_rank == 0) then                        ! This is a _big_ IF
            ! THRESH
            write (stderr_unit,*)"An oscillation in the large-component of the ",&
                "radial wavefunction is diregarded"
            write (stderr_unit, *) "for the purposes of node counting if its ", &
                "amplitude is less than 1/20 the"
            write (stderr_unit, *) "maximum amplitude.   Revise this?"
            yes = read_yes_no()
            if (yes) then
                do
                    write (stderr_unit, *) "Enter the new threshold value:"
                    read (*, *) relative_extremum_cutoff
                    if (relative_extremum_cutoff <= 0.0_real64) then
                        write (stderr_unit, *) routine_name, ": This must exceed 0;"
                    else
                        exit
                    end if
                end do
            end if
            yes = .false.
            ! METHOD

            ! Piece only for printing...

            do i = 1, 4

                do j = 1, orbital_count
                    if (.not.(radial_solution_methods(j)==i .and. .not.is_orbital_fixed(j))) cycle
                    write (stderr_unit, *) "Method ", i, " is used for ", &
                        "integrating the radial differential ", &
                        "equation for subshells"
                    exit
                end do
                if (j > orbital_count) cycle

                iend = 0
                do j = 1, orbital_count
                    if (radial_solution_methods(j)==i .and. .not.is_orbital_fixed(j)) then
                        ibeg = iend + 1
                        iend = ibeg
                        record(ibeg:iend) = " "
                        call format_integer (orbital_principal_n(j), number_text, text_length)
                        ibeg = iend + 1
                        iend = ibeg + text_length - 1
                        record(ibeg:iend) = number_text(1:text_length)
                        ibeg = iend + 1
                        if (orbital_kappa(j) < 0) then
                            iend = ibeg
                            record(ibeg:iend) = orbital_symmetry_labels(j)(1:1)
                        else
                            iend = ibeg + 1
                            record(ibeg:iend) = orbital_symmetry_labels(j)(1:2)
                        end if
                    end if
                    if (iend <= 76) cycle
                    write (stderr_unit, *) record(1:iend)
                    iend = 0
                end do
                if (iend<=0 .or. mpi_rank/=0) cycle
                write (stderr_unit, *) record(1:iend)
            end do

            ! Reads user inputs and fills array indx(1:nsubs) where nsubs itself
            ! is an output from getrsl.
            ! METHOD(1:4) is the only output. indx() and nsubs are discarded

            write (stderr_unit, *) "Select a different integration method for ", &
                "any subshell radial wavefunction?"
            yes = read_yes_no()
            if (yes) then
                do i = 1, 4
                    write (stderr_unit, *) "Method ", i, ":"
                    call select_subshell_labels (indx, nsubs)
                    do j = 1, nsubs
                        loc = indx(j)
                        if (is_orbital_fixed(loc)) cycle
                        radial_solution_methods(loc) = i
                    end do
                end do
            end if
            ! NOINVT
            write (stderr_unit, *) "The first oscillation of the large component"

            if (any(suppress_radial_inversion(:orbital_count) .and. .not.is_orbital_fixed(:orbital_count))) then
                write (stderr_unit, *) "of the following radial wavefunctions ", &
                    "will be required to be positive"
                iend = 0
                do i = 1, orbital_count
                    if (suppress_radial_inversion(i) .and. .not.is_orbital_fixed(i)) then
                        ibeg = iend + 1
                        iend = ibeg
                        record(ibeg:iend) = " "
                        call format_integer (orbital_principal_n(i), number_text, text_length)
                        ibeg = iend + 1
                        iend = ibeg + text_length - 1
                        record(ibeg:iend) = number_text(1:text_length)
                        ibeg = iend + 1
                        if (orbital_kappa(i) < 0) then
                            iend = ibeg
                            record(ibeg:iend) = orbital_symmetry_labels(i)(1:1)
                        else
                            iend = ibeg + 1
                            record(ibeg:iend) = orbital_symmetry_labels(i)(1:2)
                        end if
                    end if
                    if (iend <= 76) cycle
                    write (stderr_unit, *) record(1:iend)
                    iend = 0
                end do
                if (iend > 0) write (stderr_unit, *) record(1:iend)
                write (stderr_unit, *) "Revise this?"
            else
                write (stderr_unit, *) "of all radial wavefunctions will be required ", &
                    "to be positive.   Revise this?"
            end if
            yes = read_yes_no()
            if (yes) then
                write (stderr_unit, *) "Suppressing enforcement of positive first ", &
                    "oscillation:"
                call select_subshell_labels (indx, nsubs)
                do i = 1, nsubs
                    loc = indx(i)
                    if (is_orbital_fixed(loc)) cycle
                    suppress_radial_inversion(loc) = .true.
                end do
            end if
            ! ODAMP
            if (any(orbital_damping_factors(:orbital_count)/=0.0_real64 .and. .not.is_orbital_fixed(:orbital_count))) then
                write (stderr_unit, *) "Subshell accelerating parameters have ", &
                    "been set.   Revise these?"
            else
                write (stderr_unit, *) "Set accelerating parameters for subshell ", &
                    "radial wavefunctions?"
            end if
            yes = read_yes_no()
            if (yes) then
                write (stderr_unit, *) "Different accelerating parameters for ", &
                    "different subshell radial wavefunction?"
                yes = read_yes_no()
                if (yes) then
                    do
                        write (stderr_unit, *) "Enter an accelerating parameter"
                        write (stderr_unit, *) " (0< ODAMP < 1 allows ODAMP to be ", &
                            "reduced as convergence is approached;"
                        write (stderr_unit, *) " -1 < ODAMP < 0 implies |ODAMP| is ", &
                            "held constant):"
                        read (*, *) odampu
                        if (abs(odampu)==0.0_real64 .or. abs(odampu)>=1.0_real64) then
                            write (stderr_unit, *) routine_name, ": Value out of range ..."
                            cycle
                        end if
                        call select_subshell_labels (indx, nsubs)
                        do i = 1, nsubs
                            loc = indx(i)
                            if (is_orbital_fixed(loc)) cycle
                            orbital_damping_factors(loc) = odampu
                        end do
                        exit
                    end do
                else
                    do
                        write (stderr_unit, *) "Enter the accelerating parameter"
                        write (stderr_unit, *) " (0< ODAMP < 1 allows ODAMP to be ", &
                            "reduced as convergence is approached;"
                        write (stderr_unit, *) " -1 < ODAMP < 0 implies |ODAMP| is ", &
                            "held constant):"
                        read (*, *) odampu
                        if (abs(odampu)==0.0_real64 .or. abs(odampu)>=1.0_real64) then
                            write (stderr_unit, *) routine_name, ": Value out of range ..."
                            cycle
                        end if
                        where (.not.is_orbital_fixed(:orbital_count))
                            orbital_damping_factors(:orbital_count) = odampu
                        end where
                        exit
                    end do
                end if
            end if
            ! CDAMP
            write (stderr_unit, *) "Set accelerating parameters for the ", &
                "eigenvectors?"
            yes = read_yes_no()
            if (yes) then
                write (stderr_unit, *) "Different accelerating parameters ", &
                    "for each eigenvector?"
                yes = read_yes_no()
                if (yes) then
                    write (stderr_unit, *) "Enter an accelerating parameter for"
                    call format_integer (selected_level_count, record, text_length)
                    write (stderr_unit, *) " each of the "//record(1:text_length)//" levels :"
                    read (*, *) (mixing_damping_factors(i),i=1,selected_level_count)
                else
                    write (stderr_unit, *) "Enter the accelerating parameter:"
                    read (*, *) cdampu
                    mixing_damping_factors(:selected_level_count) = cdampu
                end if
            end if
            ! NSIC
            write (stderr_unit, *) "Following the improvement of each of the ", &
                "subshell radial wavefunctions in turn, "
            write (stderr_unit, *) "the ", orbital_updates_per_scf_iteration, " least self-consistent", &
                " functions will be improved at the"
            write (stderr_unit, *) "end of the first SCF cycle. Revise this ", &
                "setting?"
            yes = read_yes_no()
            if (yes) then
                write (stderr_unit, *) "Enter the number of additional ", "improvements:"
                read (*, *) orbital_updates_per_scf_iteration
            end if
            ! NSOLV
            write (stderr_unit, *) "The maximum number of cycles in attempting ", &
                "to solve each radial equation is "
            write (stderr_unit, *) minimum_radial_solver_iterations, " times the principal quantum", &
                " number of the radial"
            write (stderr_unit, *) "wave-function to be estimated.   ", &
                "Revise this setting?"
            yes = read_yes_no()
            if (yes) then
                write (stderr_unit, *) "Enter the factor that multiplies the ", &
                    "principal quantum number:"
                read (*, *) minimum_radial_solver_iterations
            end if
            ! Orthogonalisation
            if (enforce_orbital_orthogonality) then
                write (stderr_unit, *) "Subshell radial wavefunctions will be ", &
                    "Schmidt orthogonalised immediately"
                write (stderr_unit, *) "following their estimation to all ", &
                    "functions with poorer self-consistency."
                write (stderr_unit, *) " Revise this?"
                yes = read_yes_no()
                if (yes) enforce_orbital_orthogonality = .false.
            else
                write (stderr_unit, *) "Subshell radial wavefunctions will be ", &
                    "Schmidt orthogonalised at the end of"
                write (stderr_unit, *) "each SCF cycle.   Revise this?"
                yes = read_yes_no()
                if (yes) enforce_orbital_orthogonality = .true.
            end if
        end if                                      ! end of the _big_ IF
        return
    end subroutine read_scf_options
end module read_scf_options_interface
