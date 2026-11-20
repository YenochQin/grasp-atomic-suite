!***********************************************************************
!                                                                      *
program ris4_mpi_main
    ! Module dependencies
    !                                                                      *
    !   Entry routine for the MPI RIS program.                             *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use calculation_input_defaults, only: nondefault_options_flag
    use standard_io_units, only: stderr_unit
    use debug_output_settings,         only: cutoff
    use mpi_state, only: mpi_bcast, mpi_character, mpi_comm_world, mpi_error_code, mpi_finalize, mpi_integer, &
        mpi_process_count, mpi_rank, processor_name, processor_name_length
    use calculation_state_files, only: read_state_name, property_state_arguments
    use yes_no_input, only: read_yes_no

    ! Procedure interfaces
    use configure_debug_output_interface, only: configure_debug_output

    ! Module dependencies
    use machine_constants, only: initialize_machine_constants
    use physical_constants, only: initialize_physical_constants

    ! Procedure interfaces
    use open_summary_file_interface, only: open_summary_file

    ! Module dependencies
    use mpi_csf_basis_input, only: load_csf_basis_mpi

    ! Procedure interfaces
    use read_mixing_blocks_interface, only: read_mixing_blocks
    use prepare_isotope_calculation_interface, only: prepare_isotope_calculation

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table

    ! Procedure interfaces
    use calculate_isotope_shifts_mpi_interface, only: calculate_isotope_shifts_mpi

    ! Module dependencies
    use mpi_utilities, only: mpix_shutdown
    use mpi_utilities, only: mpix_startup
    implicit none (type, external)
    integer :: input_status, mixing_file_kind, unused_core_subshell_count, start_clock_count
    logical :: yes, use_args
    character :: name*24

    call mpix_startup(mpi_rank, mpi_process_count, processor_name, processor_name_length, start_clock_count, &
        "RIS4_MPI", "This is the MPI RIS program", &
        "isodata, descriptor Parquet/header, name.mix/name.cimix, name.wfn", &
        "name.i or name.ci", .true.)

    name = " "
    mixing_file_kind = 1
    nondefault_options_flag = 0
    use_args = .false.
    cutoff = 1.0e-10_real64

    if (mpi_rank == 0) then
        use_args = command_argument_count() > 0
        if (use_args) then
            call property_state_arguments("ris4_mpi", name, mixing_file_kind, input_status)
            if (input_status /= 0) nondefault_options_flag = -input_status
        end if
        if (.not.use_args .and. nondefault_options_flag == 0) then
            write (stderr_unit, *)
            write (stderr_unit, *) "Default settings?"
            yes = read_yes_no()
            write (stderr_unit, *)
            if (yes) then
                nondefault_options_flag = 0
            else
                nondefault_options_flag = 1
                write (stderr_unit, *) "RIS4_MPI currently supports default settings only."
            end if
        end if

        if (nondefault_options_flag == 0 .and. .not.use_args) then
            call read_state_name(name, status=input_status)
            if (input_status /= 0) then
                nondefault_options_flag = -input_status
            else
                write (stderr_unit, *)
                write (stderr_unit, *) "Mixing coefficients from a CI calc.?"
                yes = read_yes_no ()
                if (yes) then
                    mixing_file_kind = 0
                else
                    mixing_file_kind = 1
                end if
                write (stderr_unit, *)
            end if
        end if
    end if

    call mpi_bcast(nondefault_options_flag, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
    if (nondefault_options_flag /= 0) then
        call mpi_finalize(mpi_error_code)
        if (nondefault_options_flag == -2) stop 2
        stop
    end if
    call mpi_bcast(name, len(name), mpi_character, 0, mpi_comm_world, mpi_error_code)
    call mpi_bcast(mixing_file_kind, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)

    call configure_debug_output
    call initialize_machine_constants
    call initialize_physical_constants
    if (mpi_rank == 0) call open_summary_file(name,mixing_file_kind)
    call load_csf_basis_mpi(name,unused_core_subshell_count)
    call prepare_isotope_calculation(name)
    call read_mixing_blocks(name,mixing_file_kind)
    call initialize_log_factorial_table
    call calculate_isotope_shifts_mpi(name)

    call mpix_shutdown(mpi_rank, start_clock_count, "RIS4_MPI")
    stop
end program ris4_mpi_main
