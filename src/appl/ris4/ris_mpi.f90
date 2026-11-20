!***********************************************************************
!                                                                      *
      program ris_mpi
!                                                                      *
!   Entry routine for the MPI RIS program.                             *
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64
      use calculation_input_defaults, only: inpci, ndef, ndump
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use debug_output_settings,         only: cutoff
      use mpi_state, only: mpi_bcast, mpi_character, mpi_comm_world, mpi_error_code, mpi_finalize, mpi_integer, &
         & mpi_process_count, mpi_rank, processor_name, processor_name_length
      use calculation_state_files, only: read_state_name, property_state_arguments
      use yes_no_input, only: getyn
      use setdbg_i, only: setdbg
      use machine_constants, only: setmc
      use physical_constants, only: setcon
      use setsum_i, only: setsum
      use mpi_csf_basis_input, only: setcsla => load_csf_basis_mpi
      use getmixblock_i, only: getmixblock
      use getsmd_i, only: getsmd
      use log_factorial_table, only: factt
      use ris_cal_mpi_i, only: ris_cal_mpi
      use mpi_utilities, only: mpix_shutdown
      use mpi_utilities, only: mpix_startup
      implicit none (type, external)
      integer :: input_status, nci, ncore_not_used, ncount1
      logical :: yes, use_args
      character :: name*24

      call mpix_startup(mpi_rank, mpi_process_count, processor_name, processor_name_length, ncount1, &
         "RIS4_MPI", "This is the MPI RIS program", &
         "isodata, descriptor Parquet/header, name.mix/name.cimix, name.wfn", &
         "name.i or name.ci", .true.)

      name = " "
      nci = 1
      ndef = 0
      use_args = .false.
      cutoff = 1.0e-10_real64

      if (mpi_rank == 0) then
         use_args = command_argument_count() > 0
         if (use_args) then
            call property_state_arguments("ris4_mpi", name, nci, input_status)
            if (input_status /= 0) ndef = -input_status
         end if
         if (.not.use_args .and. ndef == 0) then
            write (stderr_unit, *)
            write (stderr_unit, *) "Default settings?"
            yes = getyn()
            write (stderr_unit, *)
            if (yes) then
               ndef = 0
            else
               ndef = 1
               write (stderr_unit, *) "RIS4_MPI currently supports default settings only."
            end if
         end if

         if (ndef == 0 .and. .not.use_args) then
            call read_state_name(name, status=input_status)
            if (input_status /= 0) then
               ndef = -input_status
            else
               write (stderr_unit, *)
               write (stderr_unit, *) "Mixing coefficients from a CI calc.?"
               yes = getyn ()
               if (yes) then
                  nci = 0
               else
                  nci = 1
               end if
               write (stderr_unit, *)
            end if
         end if
      end if

      call mpi_bcast(ndef, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
      if (ndef /= 0) then
         call mpi_finalize(mpi_error_code)
         if (ndef == -2) stop 2
         stop
      end if
      call mpi_bcast(name, len(name), mpi_character, 0, mpi_comm_world, mpi_error_code)
      call mpi_bcast(nci, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)

      call setdbg
      call setmc
      call setcon
      if (mpi_rank == 0) call setsum(name,nci)
      call setcsla(name,ncore_not_used)
      call getsmd(name)
      call getmixblock(name,nci)
      call factt
      call ris_cal_mpi(name)

      call mpix_shutdown(mpi_rank, ncount1, "RIS4_MPI")
      stop
      end program ris_mpi
