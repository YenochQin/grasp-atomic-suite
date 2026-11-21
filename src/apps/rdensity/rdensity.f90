! Shared entry point for rdensity and rdensity_mpi.
program rdensity
    ! Module dependencies
    use calculation_state_files, only: orbital_input_path, mixing_input_path, read_state_name, property_state_arguments
    use calculation_input_defaults, only: nondefault_options_flag
    use debug_output_settings, only: matrix_element_screening_cutoff
    use density_parallel, only: broadcast_options, density_fail, parallel_finish, parallel_start, mpi_rank, mpi_process_count
    use density_calculation, only: calculate_density
    use yes_no_input, only: read_yes_no

    ! Procedure interfaces
    use configure_debug_output_interface, only: configure_debug_output

    ! Module dependencies
    use machine_constants, only: initialize_machine_constants
    use physical_constants, only: initialize_physical_constants
    use density_parallel, only: load_density_basis

    ! Procedure interfaces
    use read_radial_density_input_interface, only: read_radial_density_input
    use read_mixing_blocks_interface, only: read_mixing_blocks

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    character(24) :: name
    character(256) :: command_argument, io_message
    integer :: mixing_file_kind, ordering, status, io_status, core_subshell_count
    logical :: yes, exists, interactive

    call parallel_start
    name = ""
    mixing_file_kind = 1
    nondefault_options_flag = 0
    ordering = 2
    status = 0
    matrix_element_screening_cutoff = 1.0e-10_real64
    interactive = command_argument_count() == 0
    if (mpi_rank == 0) then
        if (interactive) then
            write(*,*) "Default settings?"
            yes = read_yes_no()
            if (.not.yes) nondefault_options_flag = 1
            if (mpi_process_count > 1 .and. nondefault_options_flag /= 0) then
                write(*,*) "MPI runs support default settings only, as in RHFS/RIS."
                status = 2
            else
                call read_state_name(name, status=status)
                if (status == 0) then
                    write(*,*) "Mixing coefficients from a CI calc.?"
                    yes = read_yes_no()
                    if (yes) mixing_file_kind = 0
                end if
            end if
        else
            call property_state_arguments("rdensity[_mpi]", name, mixing_file_kind, status, ordering)
        end if
    end if
    call broadcast_options(name, mixing_file_kind, nondefault_options_flag, ordering, status)
    if (status == 1) then
        call parallel_finish
        stop
    end if
    if (status /= 0) call density_fail("Invalid command line or interactive input")

    call require_file("isodata")
    call require_file(orbital_input_path(name))
    call require_file(mixing_input_path(name, mixing_file_kind == 0))
    if (mixing_file_kind == 0) then
        command_argument = trim(name)//".cd"
    else
        command_argument = trim(name)//".d"
    end if
    if (mpi_rank == 0) then
        inquire(file=trim(command_argument), exist=exists)
        if (exists) call density_fail("Density output already exists: "//trim(command_argument))
    end if
    call configure_debug_output
    call initialize_machine_constants
    call initialize_physical_constants
    ! Read-only inputs are shared; angular and orbital state is rank-local.
    call load_density_basis(name, core_subshell_count)
    call read_radial_density_input(name)
    call read_mixing_blocks(name, mixing_file_kind)
    call initialize_log_factorial_table
    if (mpi_rank == 0 .and. interactive) then
        write(*,*) "How do you want to order your eigenvectors?"
        write(*,*) "1) Dominant component; 2) Decreasing occupation (preserving pure orbitals)"
        read(*,*,iostat=io_status,iomsg=io_message) ordering
        if (io_status /= 0) then
            write(*,*) "Cannot read orbital ordering: ", trim(io_message)
            status = 2
        else if (ordering /= 1 .and. ordering /= 2) then
            status = 2
        end if
    end if
    call broadcast_options(name, mixing_file_kind, nondefault_options_flag, ordering, status)
    if (status /= 0) call density_fail("Orbital ordering must be 1 or 2")
    call calculate_density(name, mixing_file_kind, ordering)
    call parallel_finish
contains
    subroutine require_file(filename)
        ! Arguments
        character(*), intent(in) :: filename

        ! Local declarations
        logical :: found
        inquire(file=filename, exist=found)
        if (.not.found) call density_fail("Missing input file: "//filename)
    end subroutine require_file
end program rdensity
