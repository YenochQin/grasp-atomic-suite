!***********************************************************************
!                                                                      *
      program gj_mpi
!                                                                      *
!   Entry routine for the MPI Landé-factor program.                    *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use default_c, only: double, inpci, ndef, ndump
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      use memory_man, only: alloc, dalloc, double
      use mpi_c, only: double, mpi_bcast, mpi_character, mpi_comm_world, mpi_error_code, mpi_finalize, mpi_integer, &
         & mpi_process_count, mpi_rank, processor_name, processor_name_length
      use prnt_c,          only: nvec
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use compute_gj_only_mpi_i, only: compute_gj_only_mpi
      use gethfd_i, only: gethfd
      use getmixblock_i, only: getmixblock
      use getyn_i, only: getyn
      use setcon_i, only: setcon
      use csf_basis_mpi, only: setcsla => load_csf_basis_mpi
      use setdbg_i, only: setdbg
      use setmc_i, only: setmc
      use setout_gj_i, only: setout_gj
      use factt_i, only: factt
      use mpiu_i, only: mpix_shutdown
      use mpiu_i, only: mpix_startup
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: argc, iarg, k, nci, ncore_not_used, ncount1
      logical :: yes, use_args
      real(real64), pointer :: gjc_diag(:), dgjc_diag(:)
      character(len=24) :: arg
      character :: name*24
!-----------------------------------------------
!
      call mpix_startup(mpi_rank, mpi_process_count, processor_name, processor_name_length, ncount1, &
         "GJ_MPI", "This is the MPI Lande-factor program", &
         "isodata, descriptor Parquet/header, name.(c)m, name.w", &
         "name.gj or name.cgj", .true.)
!
      name = " "
      nci = 1
      ndef = 0
      use_args = .false.
!
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
                  write (stderr_unit, *) "Usage: gj_mpi [name] [--ci|--nonci]"
                  write (stderr_unit, *) "If no arguments are given, interactive mode is used."
                  ndef = -1
               else if (iarg /= 1 .or. arg(1:2) == "--") then
                  write (stderr_unit, *) "Unrecognized argument: ", trim(arg)
                  write (stderr_unit, *) "Usage: gj_mpi [name] [--ci|--nonci]"
                  ndef = -1
               end if
            end do
            if (ndef == 0 .and. len_trim(name) == 0) then
               write (stderr_unit, *) "Usage: gj_mpi [name] [--ci|--nonci]"
               ndef = -1
            end if
         end if
!
         if (.not.use_args .and. ndef == 0) then
            write (stderr_unit, *)
            write (stderr_unit, *) "Default settings?"
            yes = getyn()
            write (stderr_unit, *)
            if (yes) then
               ndef = 0
            else
               ndef = 1
               write (stderr_unit, *) "GJ_MPI currently supports default settings only."
            end if
         end if
!
         if (ndef == 0 .and. .not.use_args) then
   10       continue
            write (stderr_unit, *) "Name of state"
            read (*, "(A)") name
            k = index(name," ")
            if (k == 1) then
               write (stderr_unit, *) "Names may not start with a blank"
               goto 10
            end if
!
            write (stderr_unit, *)
            write (stderr_unit, *) "Mixing coefficients from a CI calc.?"
            yes = getyn()
            if (yes) then
               nci = 0
            else
               nci = 1
            end if
         end if
      end if
!
      call mpi_bcast(ndef, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
      if (ndef /= 0) then
         call mpi_finalize(mpi_error_code)
         stop
      end if
!
      call mpi_bcast(name, len(name), mpi_character, 0, mpi_comm_world, mpi_error_code)
      call mpi_bcast(nci, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
!
      call setdbg
      call setmc
      call setcon
      call setcsla(name, ncore_not_used)
      call gethfd(name)
      call getmixblock(name, nci)
      call factt
!
      call alloc(gjc_diag, nvec, "GJC_DIAG", "GJ_MPI")
      call alloc(dgjc_diag, nvec, "DGJC_DIAG", "GJ_MPI")
!
      call compute_gj_only_mpi(gjc_diag, dgjc_diag)
      if (mpi_rank == 0) call setout_gj(name, nci, gjc_diag, dgjc_diag)
!
      call dalloc(gjc_diag, "GJC_DIAG", "GJ_MPI")
      call dalloc(dgjc_diag, "DGJC_DIAG", "GJ_MPI")
!
      call mpix_shutdown(mpi_rank, ncount1, "GJ_MPI")
      stop
      end program gj_mpi
