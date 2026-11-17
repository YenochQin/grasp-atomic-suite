!***********************************************************************
!                                                                      *
      program ris_mpi
!                                                                      *
!   Entry routine for the MPI RIS program.                             *
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64
      use default_c, only: double, inpci, ndef, ndump
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      use debug_c,         only: cutoff
      use mpi_c, only: double, mpi_bcast, mpi_character, mpi_comm_world, mpi_error_code, mpi_finalize, mpi_integer, &
         & mpi_process_count, mpi_rank, processor_name, processor_name_length
      use getyn_i, only: getyn
      use setdbg_i, only: setdbg
      use setmc_i, only: setmc
      use setcon_i, only: setcon
      use setsum_i, only: setsum
      use csf_basis_mpi, only: setcsla => load_csf_basis_mpi
      use getmixblock_i, only: getmixblock
      use getsmd_i, only: getsmd
      use factt_i, only: factt
      use ris_cal_mpi_i, only: ris_cal_mpi
      use mpiu_i, only: mpix_shutdown
      use mpiu_i, only: mpix_startup
      implicit none (type, external)
      integer :: argc, iarg, k, nci, ncore_not_used, ncount1
      logical :: yes, use_args
      character(len=24) :: arg
      character :: name*24

      call mpix_startup(mpi_rank, mpi_process_count, processor_name, processor_name_length, ncount1, &
         "RIS4_MPI", "This is the MPI RIS program", &
         "isodata, descriptor Parquet/header, name.(c)m, name.w", &
         "name.i or name.ci", .true.)

      name = " "
      nci = 1
      ndef = 0
      use_args = .false.
      cutoff = 1.0e-10_real64

      if (mpi_rank == 0) then
         argc = command_argument_count()
         if (argc > 0) then
            use_args = .true.
            do iarg = 1, argc
               arg = " "
               call get_command_argument(iarg, arg)
               if (iarg == 1 .and. arg(1:2) /= "--") then
                  name = arg
               else if (trim(arg) == "--ci") then
                  nci = 0
               else if (trim(arg) == "--nonci") then
                  nci = 1
               else if (trim(arg) == "--help" .or. trim(arg) == "-h") then
                  write (stderr_unit, *) "Usage: ris4_mpi [name] [--ci|--nonci]"
                  write (stderr_unit, *) "If no arguments are given, interactive mode is used."
                  ndef = -1
               else
                  write (stderr_unit, *) "Unrecognized argument: ", trim(arg)
                  write (stderr_unit, *) "Usage: ris4_mpi [name] [--ci|--nonci]"
                  ndef = -1
               end if
            end do
            if (ndef == 0 .and. len_trim(name) == 0) then
               write (stderr_unit, *) "Usage: ris4_mpi [name] [--ci|--nonci]"
               ndef = -1
            end if
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
    9       write (stderr_unit, *) "Name of state"
            read(*,"(A)") name
            k=index(name," ")
            if (k==1) then
               write (stderr_unit, *) "Names may not start with a blank"
               goto 9
            end if
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

      call mpi_bcast(ndef, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
      if (ndef /= 0) then
         call mpi_finalize(mpi_error_code)
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
