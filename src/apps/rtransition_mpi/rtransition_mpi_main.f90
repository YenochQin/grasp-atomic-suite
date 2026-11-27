! Calculate radiative multipole transitions using MPI-distributed angular contributions.
! Use orbitals transformed by rbiotransform so the two orbital sets are biorthogonal.
! Written by Per Jonsson,   Department of Computer Science
!                         Vanderbilt University, USA
!
! Modified by Gediminas Gaigalas for new spin-angular integration
! and for reducing usage of CPU memory.        NIST, October 2017
program rtransition_mpi_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    use calculation_input_defaults, only: input_mixing_file_kind, nondefault_options_flag, angular_output_dump_flag
    use yes_no_input, only: read_yes_no
    use machine_constants, only: initialize_machine_constants
    use physical_constants, only: initialize_physical_constants

    ! Procedure interfaces
    use mpi_calculation_state, only: read_mpi_state_pair
    use merge_transition_csf_bases_interface, only: merge_transition_csf_bases

    use read_transition_options_interface, only: read_transition_options
    use write_input_summary_interface, only: write_input_summary

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table

    ! Procedure interfaces
    use calculate_oscillator_strengths_interface, only: calculate_oscillator_strengths

    ! Module dependencies
    use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_process_count, mpi_rank, &
        processor_name, processor_name_length

    use posix_services, only: terminate_process
    use mpi_utilities, only: initialize_mpi_runtime_with_status
    use mpi_utilities, only: finish_mpi_application
    use isotope_data_input, only: read_isotope_data
    use mpi_change_directory, only: change_mpi_directory
    implicit none (type, external)
    integer :: ntest,label_output_option,permanent_directory_length,scratch_directory_length, iii
    integer(int64) :: start_clock_count
    logical :: yes
    character :: name(2)*24
    character :: fullname(2)*128
    character(len=3)  ::  block_label_text
    character(len=128) :: isotope_file_path,startup_directory,permanent_directory,scratch_directory
    ! Initialize MPI and set up the I/O and scratch directories.
    call initialize_mpi_runtime_with_status (mpi_rank, mpi_process_count, processor_name, processor_name_length, &
        & start_clock_count,           &
        startup_directory, permanent_directory, scratch_directory, "RTRANSITION_MPI")
    write (block_label_text, "(I3.3)") mpi_rank
    print*, scratch_directory , " = tmpdir"

    permanent_directory_length = len_trim (permanent_directory)
    scratch_directory_length = len_trim (scratch_directory)

    isotope_file_path = trim(startup_directory)//"/isodata"
    ntest = 1001

    if (mpi_rank == 0) then
        write (output_unit, *)
        write (output_unit, *) " Default settings?"
        yes = read_yes_no()
        write (output_unit, *)
        if (yes) then
            nondefault_options_flag = 0
            angular_output_dump_flag = 1
        else
            nondefault_options_flag = 1
            write (output_unit, *) " Dump angular data to file?"
            yes = read_yes_no()
            if (yes) then
                angular_output_dump_flag = 1
            else
                angular_output_dump_flag = 0
            end if
        end if
        write (output_unit, *)
        write (output_unit, *) " Input from a CI calculation?"
        yes = read_yes_no()
        write (output_unit, *)
        if (yes) then
            input_mixing_file_kind = 0
        else
            input_mixing_file_kind = 1
        end if
    end if   !myid=0
    call mpi_bcast (nondefault_options_flag,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
    call mpi_bcast (angular_output_dump_flag,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
    call mpi_bcast (input_mixing_file_kind,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
    ! start
    ! end
    !
    ! Initialize machine-dependent constants.
    call initialize_machine_constants
    ! Initialize physical constants.
    call initialize_physical_constants
    ! Obtain the names of the initial and final state files
    call read_mpi_state_pair(name)
    fullname(1) = trim(startup_directory) // "/" // trim(name(1))
    fullname(2) = trim(startup_directory) // "/" // trim(name(2))
    ! Load and align the initial and final CSF bases in memory.
    if (mpi_rank == 0) then
        print*,"FULLNAME1=",fullname(1)
        print*,"FULLNAME2=",fullname(2)
    end if

    call merge_transition_csf_bases (fullname)
    !
    !
    !
    !
    ! Read the nuclear parameters from isodata.
    call read_isotope_data(isotope_file_path)

    ! Get the remaining information
    call read_transition_options (fullname)
    ! Open the transition report and record the input settings.
    if (mpi_rank == 0) then
        iii = len_trim(startup_directory)
        call change_mpi_directory(trim(startup_directory),iii,mpi_error_code)
        if (mpi_error_code/=0) then
            print *, "error changing dir!"
            call terminate_process(1)
        end if
        label_output_option = 0
        call write_input_summary(name,input_mixing_file_kind,label_output_option)
        iii = len_trim(scratch_directory)
        call change_mpi_directory(trim(scratch_directory)//"/"//block_label_text,iii+4,mpi_error_code)
        if (mpi_error_code/=0) then
            print *, "Error changing dir!"
            call terminate_process(1)
        end if
    end if  !myid = 0
    ! Initialize the log-factorial table for angular algebra.
    call initialize_log_factorial_table
    ! Proceed with the transition calculation
    call calculate_oscillator_strengths (name,fullname,scratch_directory,startup_directory,block_label_text)

    call finish_mpi_application (mpi_rank, start_clock_count, "RTRANSITION_MPI")

    stop
end program rtransition_mpi_main
