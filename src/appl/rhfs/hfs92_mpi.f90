!***********************************************************************
!                                                                      *
      program hfs92_mpi
!                                                                      *
!   Entry routine for the MPI hyperfine-structure program.             *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use default_c, only: double, inpci, ndef, ndump
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      use mpi_c, only: double, mpi_bcast, mpi_character, mpi_comm_world, mpi_error_code, mpi_finalize, mpi_integer, &
         & mpi_process_count, mpi_rank, processor_name, processor_name_length
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i, only: getyn
      use setdbg_i, only: setdbg
      use setmc_i, only: setmc
      use setcon_i, only: setcon
      use setsum_i, only: setsum
      use setcsla_i, only: setcsla
      use gethfd_i, only: gethfd
      use getmixblock_i, only: getmixblock
      use strsum_i, only: strsum
      use factt_i, only: factt
      use hfsgg_mpi_i, only: hfsgg_mpi
      use mpiu_i, only: gdmaxmpi_root
      use mpiu_i, only: mpix_shutdown
      use mpiu_i, only: mpix_startup
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: argc, iarg, k, nci, ncore_not_used, ncount1
      integer :: nstage1, nstage2, ncount_rate, ncount_max
      logical :: yes, use_args
      real(real64) :: time_setcsla, time_gethfd, time_getmixblock, time_hfsgg
      real(real64) :: timbuf(1)
      character(len=24) :: arg
      character :: name*24
!-----------------------------------------------
!
      call mpix_startup(mpi_rank, mpi_process_count, processor_name, processor_name_length, ncount1, &
         "RHFS_MPI", "This is the MPI hyperfine structure program", &
         "isodata, name.c, name.(c)m, name.w", &
         "name.(c)h, name.(c)hoffd", .true.)
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
                  write (stderr_unit, *) "Usage: rhfs_mpi [name] [--ci|--nonci]"
                  write (stderr_unit, *) "If no arguments are given, interactive mode is used."
                  ndef = -1
               else if (iarg /= 1 .or. arg(1:2) == "--") then
                  write (stderr_unit, *) "Unrecognized argument: ", trim(arg)
                  write (stderr_unit, *) "Usage: rhfs_mpi [name] [--ci|--nonci]"
                  ndef = -1
               end if
            end do
            if (ndef == 0 .and. len_trim(name) == 0) then
               write (stderr_unit, *) "Usage: rhfs_mpi [name] [--ci|--nonci]"
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
               write (stderr_unit, *) "RHFS_MPI currently supports default settings only."
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
      if (mpi_rank == 0) call setsum(name, nci)
      call system_clock(nstage1, ncount_rate, ncount_max)
      call setcsla(name, ncore_not_used)
      call system_clock(nstage2, ncount_rate, ncount_max)
      time_setcsla = dble(nstage2 - nstage1) / dble(ncount_rate)
      timbuf(1) = time_setcsla
      call gdmaxmpi_root(timbuf, 1)
      if (mpi_rank == 0) time_setcsla = timbuf(1)
!
      call system_clock(nstage1, ncount_rate, ncount_max)
      call gethfd(name)
      call system_clock(nstage2, ncount_rate, ncount_max)
      time_gethfd = dble(nstage2 - nstage1) / dble(ncount_rate)
      timbuf(1) = time_gethfd
      call gdmaxmpi_root(timbuf, 1)
      if (mpi_rank == 0) time_gethfd = timbuf(1)
!
      call system_clock(nstage1, ncount_rate, ncount_max)
      call getmixblock(name, nci)
      call system_clock(nstage2, ncount_rate, ncount_max)
      time_getmixblock = dble(nstage2 - nstage1) / dble(ncount_rate)
      timbuf(1) = time_getmixblock
      call gdmaxmpi_root(timbuf, 1)
      if (mpi_rank == 0) time_getmixblock = timbuf(1)
      if (mpi_rank == 0) call strsum
      call factt
!
      call system_clock(nstage1, ncount_rate, ncount_max)
      call hfsgg_mpi
      call system_clock(nstage2, ncount_rate, ncount_max)
      time_hfsgg = dble(nstage2 - nstage1) / dble(ncount_rate)
      timbuf(1) = time_hfsgg
      call gdmaxmpi_root(timbuf, 1)
      if (mpi_rank == 0) time_hfsgg = timbuf(1)
!
      if (mpi_rank == 0) then
         write (6, "(A, F10.3, A)") "RHFS_MPI setcsla wall time  (max rank): ", &
            time_setcsla, " s"
         write (6, "(A, F10.3, A)") "RHFS_MPI gethfd wall time   (max rank): ", &
            time_gethfd, " s"
         write (6, "(A, F10.3, A)") "RHFS_MPI getmix wall time   (max rank): ", &
            time_getmixblock, " s"
         write (6, "(A, F10.3, A)") "RHFS_MPI hfsgg wall time    (max rank): ", &
            time_hfsgg, " s"
      end if
!
      call mpix_shutdown(mpi_rank, ncount1, "RHFS_MPI")
      stop
      end program hfs92_mpi
