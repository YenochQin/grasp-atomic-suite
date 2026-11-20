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
      use calculation_input_defaults, only: inpci, ndef, ndump
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use mpi_state, only: mpi_bcast, mpi_character, mpi_comm_world, mpi_error_code, mpi_finalize, mpi_integer, &
         & mpi_process_count, mpi_rank, processor_name, processor_name_length
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use calculation_state_files, only: read_state_name, property_state_arguments
      use yes_no_input, only: getyn
      use setdbg_i, only: setdbg
      use machine_constants, only: setmc
      use physical_constants, only: setcon
      use setsum_i, only: setsum
      use mpi_csf_basis_input, only: setcsla => load_csf_basis_mpi
      use gethfd_i, only: gethfd
      use getmixblock_i, only: getmixblock
      use strsum_i, only: strsum
      use log_factorial_table, only: factt
      use hfsgg_mpi_i, only: hfsgg_mpi
      use mpi_utilities, only: gdmaxmpi_root
      use mpi_utilities, only: mpix_shutdown
      use mpi_utilities, only: mpix_startup
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: input_status, nci, unused_core_subshell_count, start_clock_count
      integer :: nstage1, nstage2, ncount_rate, ncount_max
      logical :: yes, use_args
      real(real64) :: time_setcsla, time_gethfd, time_getmixblock, time_hfsgg
      real(real64) :: timbuf(1)
      character :: name*24
!-----------------------------------------------
!
      call mpix_startup(mpi_rank, mpi_process_count, processor_name, processor_name_length, start_clock_count, &
         "RHFS_MPI", "This is the MPI hyperfine structure program", &
         "isodata, descriptor Parquet/header, name.mix/name.cimix, name.wfn", &
         "name.(c)h, name.(c)hoffd", .true.)
!
      name = " "
      nci = 1
      ndef = 0
      use_args = .false.
!
      if (mpi_rank == 0) then
         use_args = command_argument_count() > 0
         if (use_args) then
            call property_state_arguments("rhfs_mpi", name, nci, input_status)
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
               write (stderr_unit, *) "RHFS_MPI currently supports default settings only."
            end if
         end if
!
         if (ndef == 0 .and. .not.use_args) then
            call read_state_name(name, status=input_status)
            if (input_status /= 0) then
               ndef = -input_status
            else
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
      end if
!
      call mpi_bcast(ndef, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
      if (ndef /= 0) then
         call mpi_finalize(mpi_error_code)
         if (ndef == -2) stop 2
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
      call setcsla(name, unused_core_subshell_count)
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
      call mpix_shutdown(mpi_rank, start_clock_count, "RHFS_MPI")
      stop
      end program hfs92_mpi
