! Calculate hyperfine and magnetic matrix elements across MPI ranks.
program hyperfine_mpi_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use calculation_input_defaults, only: nondefault_options_flag
    use standard_io_units, only: stderr_unit
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
    use read_hyperfine_options_interface, only: read_hyperfine_options
    use read_mixing_blocks_interface, only: read_mixing_blocks
    use write_input_summary_interface, only: write_input_summary

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table

    ! Procedure interfaces
    use calculate_hyperfine_matrix_elements_mpi_interface, only: calculate_hyperfine_matrix_elements_mpi

    ! Module dependencies
    use mpi_utilities, only: reduce_mpi_real_max_to_root
    use mpi_utilities, only: finish_mpi_application
    use mpi_utilities, only: start_mpi_application
    implicit none (type, external)
    integer :: input_status, mixing_file_kind, unused_core_subshell_count
    integer(int64) :: start_clock_count
    integer :: stage_start_clock_count, stage_finish_clock_count, clock_tick_rate, maximum_clock_count
    logical :: prompt_answer, use_command_line_arguments
    real(real64) :: basis_loading_seconds, hyperfine_option_input_seconds, mixing_input_seconds, hyperfine_matrix_element_seconds
    real(real64) :: elapsed_time_reduction_buffer(1)
    character :: name*24

    call start_mpi_application(mpi_rank, mpi_process_count, processor_name, processor_name_length, start_clock_count, &
        "HYPERFINE_MPI", "This is the MPI hyperfine structure program", &
        "isodata, descriptor Parquet/header, name.mix/name.cimix, name.wfn", &
        "name.(ci)hyperfine, name.(ci)hyperfine.off-diagonal", .true.)
    name = " "
    mixing_file_kind = 1
    nondefault_options_flag = 0
    use_command_line_arguments = .false.
    if (mpi_rank == 0) then
        use_command_line_arguments = command_argument_count() > 0
        if (use_command_line_arguments) then
            call property_state_arguments("hyperfine_mpi", name, mixing_file_kind, input_status)
            if (input_status /= 0) nondefault_options_flag = -input_status
        end if
        if (.not.use_command_line_arguments .and. nondefault_options_flag == 0) then
            write (stderr_unit, *)
            write (stderr_unit, *) "Default settings?"
            prompt_answer = read_yes_no()
            write (stderr_unit, *)
            if (prompt_answer) then
                nondefault_options_flag = 0
            else
                nondefault_options_flag = 1
                write (stderr_unit, *) "HYPERFINE_MPI currently supports default settings only."
            end if
        end if
        if (nondefault_options_flag == 0 .and. .not.use_command_line_arguments) then
            call read_state_name(name, status=input_status)
            if (input_status /= 0) then
                nondefault_options_flag = -input_status
            else
                write (stderr_unit, *)
                write (stderr_unit, *) "Mixing coefficients from a CI calc.?"
                prompt_answer = read_yes_no()
                if (prompt_answer) then
                    mixing_file_kind = 0
                else
                    mixing_file_kind = 1
                end if
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
    if (mpi_rank == 0) call open_summary_file(name, mixing_file_kind)
    call system_clock(stage_start_clock_count, clock_tick_rate, maximum_clock_count)
    call load_csf_basis_mpi(name, unused_core_subshell_count)
    call system_clock(stage_finish_clock_count, clock_tick_rate, maximum_clock_count)
    basis_loading_seconds = dble(stage_finish_clock_count - stage_start_clock_count) / dble(clock_tick_rate)
    elapsed_time_reduction_buffer(1) = basis_loading_seconds
    call reduce_mpi_real_max_to_root(elapsed_time_reduction_buffer, 1)
    if (mpi_rank == 0) basis_loading_seconds = elapsed_time_reduction_buffer(1)
    call system_clock(stage_start_clock_count, clock_tick_rate, maximum_clock_count)
    call read_hyperfine_options(name)
    call system_clock(stage_finish_clock_count, clock_tick_rate, maximum_clock_count)
    hyperfine_option_input_seconds = dble(stage_finish_clock_count - stage_start_clock_count) / dble(clock_tick_rate)
    elapsed_time_reduction_buffer(1) = hyperfine_option_input_seconds
    call reduce_mpi_real_max_to_root(elapsed_time_reduction_buffer, 1)
    if (mpi_rank == 0) hyperfine_option_input_seconds = elapsed_time_reduction_buffer(1)
    call system_clock(stage_start_clock_count, clock_tick_rate, maximum_clock_count)
    call read_mixing_blocks(name, mixing_file_kind)
    call system_clock(stage_finish_clock_count, clock_tick_rate, maximum_clock_count)
    mixing_input_seconds = dble(stage_finish_clock_count - stage_start_clock_count) / dble(clock_tick_rate)
    elapsed_time_reduction_buffer(1) = mixing_input_seconds
    call reduce_mpi_real_max_to_root(elapsed_time_reduction_buffer, 1)
    if (mpi_rank == 0) mixing_input_seconds = elapsed_time_reduction_buffer(1)
    if (mpi_rank == 0) call write_input_summary
    call initialize_log_factorial_table
    call system_clock(stage_start_clock_count, clock_tick_rate, maximum_clock_count)
    call calculate_hyperfine_matrix_elements_mpi
    call system_clock(stage_finish_clock_count, clock_tick_rate, maximum_clock_count)
    hyperfine_matrix_element_seconds = dble(stage_finish_clock_count - stage_start_clock_count) / dble(clock_tick_rate)
    elapsed_time_reduction_buffer(1) = hyperfine_matrix_element_seconds
    call reduce_mpi_real_max_to_root(elapsed_time_reduction_buffer, 1)
    if (mpi_rank == 0) hyperfine_matrix_element_seconds = elapsed_time_reduction_buffer(1)
    if (mpi_rank == 0) then
        write (output_unit, "(A, F10.3, A)") "HYPERFINE_MPI setcsla wall time  (max rank): ", &
            basis_loading_seconds, " s"
        write (output_unit, "(A, F10.3, A)") "HYPERFINE_MPI gethfd wall time   (max rank): ", &
            hyperfine_option_input_seconds, " s"
        write (output_unit, "(A, F10.3, A)") "HYPERFINE_MPI getmix wall time   (max rank): ", &
            mixing_input_seconds, " s"
        write (output_unit, "(A, F10.3, A)") "HYPERFINE_MPI hfsgg wall time    (max rank): ", &
            hyperfine_matrix_element_seconds, " s"
    end if
    call finish_mpi_application(mpi_rank, start_clock_count, "HYPERFINE_MPI")
    stop
end program hyperfine_mpi_main
