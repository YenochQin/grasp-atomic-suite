!***********************************************************************
!                                                                      *
program rtransition_mpi_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    !                                                                      *
    !     This program calculates the transition parameters for a          *
    !     transition between an initial and a final state                  *
    !     The program assumes that the radial orbitals of the initial      *
    !     and final state have been transformed by the BIOTRA program      *
    !     as to become biorthonormal, in which case the normal Racah       *
    !     algebra can be used.                                             *
    !                                                                      *
    !     Written by Per Jonsson,   Department of Computer Science         *
    !                             Vanderbilt University, USA               *
    !                                                                      *
    !     Modified by Gediminas Gaigalas for new spin-angular integration  *
    !     and for reducing usage of CPU memory.        NIST, October 2017  *
    !                                                                      *
    !***********************************************************************
    use calculation_input_defaults, only: input_mixing_file_kind, nondefault_options_flag, angular_output_dump_flag
    use yes_no_input, only: read_yes_no
    use machine_constants, only: initialize_machine_constants
    use physical_constants, only: initialize_physical_constants

    ! Procedure interfaces
    use read_state_filenames_interface, only: read_state_filenames
    use merge_transition_csf_bases_interface, only: merge_transition_csf_bases

    use read_transition_options_interface, only: read_transition_options
    use write_input_summary_interface, only: write_input_summary

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table

    ! Procedure interfaces
    use calculate_oscillator_strengths_interface, only: calculate_oscillator_strengths

    ! Module dependencies
    use mpi_state, only: mpi_bcast, mpi_character, mpi_comm_world, mpi_error_code, mpi_integer, mpi_process_count, mpi_rank, &
        processor_name, processor_name_length

    use posix_services, only: terminate_process
    use mpi_utilities, only: initialize_mpi_runtime_with_status
    use mpi_utilities, only: finalize_mpi_runtime_with_status
    use isotope_data_input, only: read_isotope_data
    use mpi_change_directory, only: change_mpi_directory
    implicit none (type, external)
    integer :: ntest,start_clock_count,label_output_option,permanent_directory_length,scratch_directory_length, iii
    logical :: yes
    character :: name(2)*24
    character :: fullname(2)*128
    character(len=3)  ::  block_label_text
    character(len=128) :: isotope_file_path,startup_directory,permanent_directory,scratch_directory
    !-----------------------------------------------
    !
    !=======================================================================
    !  Start mpi --- get processor info: myid, nprocs, host name; and print
    !=======================================================================
    call initialize_mpi_runtime_with_status (mpi_rank, mpi_process_count, processor_name, processor_name_length, &
        & start_clock_count,           &
        startup_directory, permanent_directory, scratch_directory, "RTRANSITION_MPI")
    write (block_label_text, "(I3.3)") mpi_rank
    print*, scratch_directory , " = tmpdir"

    permanent_directory_length = len_trim (permanent_directory)
    scratch_directory_length = len_trim (scratch_directory)

    !=======================================================================
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
    !Rasa -- start
    !GG      LDBPR = .FALSE.
    !      WRITE (6, *) ' Generate debug output?'
    !      YES = GETYN()
    !      WRITE (6, *)
    !      IF (YES) THEN
    !         LDBPR(18) = .TRUE.
    !         WRITE (6, *) ' Enter the cutoff'
    !         READ (5, *) CUTOFF
    !      ENDIF
    !Rasa -- end
    !
    !   Perform machine- and installation-dependent setup
    !
    call initialize_machine_constants
    !
    !   Set up the physical constants
    !
    call initialize_physical_constants
    !
    !   Obtain the names of the initial and final state files
    !
    if (mpi_rank == 0) call read_state_filenames(name)
    call mpi_bcast (name,48,mpi_character,0,mpi_comm_world,mpi_error_code)
    fullname(1) = trim(startup_directory) // "/" // trim(name(1))
    fullname(2) = trim(startup_directory) // "/" // trim(name(2))
    !
    !   Open, check, load data from, and close, the initial and final state
    !   .csl  files. These files are then merged to one file.
    !
    !cjb print FULLNAME1= FULLNAME2=
    if (mpi_rank == 0) then
        print*,"FULLNAME1=",fullname(1)
        print*,"FULLNAME2=",fullname(2)
    end if

    call merge_transition_csf_bases (fullname)
    !
    !   Open, check, load data from, and close, the merged .csl  file
    !
    ! The descriptor bases were merged directly in memory by MRGCSL.
    !
    !   Read mixing coefficients
    !
    !      CALL READMIX(NAME,INPCI)
    !
    !   Test mixing coefficients
    !
    !      IF (NTEST.GT.1000) CALL TESTMIX
    !
    !   Open, check, load data from, and close the .iso file
    !
    call read_isotope_data(isotope_file_path)

    !
    !   Get the remaining information
    !
    call read_transition_options (fullname)
    !
    !   Open and append a summary of the inputs to the  .sum  file
    !
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
    !
    !   Set up the table of logarithms of factorials
    !
    call initialize_log_factorial_table
    !
    !   Proceed with the transition calculation
    !
    call calculate_oscillator_strengths (name,fullname,scratch_directory,startup_directory,block_label_text)
    !
    !   Print completion message
    !
    if (mpi_rank == 0) then
        print *
        print *, "RTRANSITION_MPI: Execution complete."
    end if
    !=======================================================================
    !  Execution finished; Statistics output
    !=======================================================================

    call finalize_mpi_runtime_with_status (mpi_rank, mpi_process_count, processor_name, processor_name_length,                     &
        start_clock_count, "RTRANSITION_MPI")
    !=======================================================================
    !
    stop
end program rtransition_mpi_main
