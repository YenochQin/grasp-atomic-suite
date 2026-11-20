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
      use calculation_input_defaults, only: inpci, ndef, ndump
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use array_memory_management, only: alloc, dalloc
      use mpi_state, only: mpi_bcast, mpi_character, mpi_comm_world, mpi_error_code, mpi_finalize, mpi_integer, &
         & mpi_process_count, mpi_rank, processor_name, processor_name_length
      use selected_level_indices,          only: nvec
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use compute_gj_only_mpi_i, only: compute_gj_only_mpi
      use gethfd_i, only: gethfd
      use getmixblock_i, only: getmixblock
      use calculation_state_files, only: read_state_name, property_state_arguments
      use yes_no_input, only: getyn
      use physical_constants, only: setcon
      use mpi_csf_basis_input, only: setcsla => load_csf_basis_mpi
      use setdbg_i, only: setdbg
      use machine_constants, only: setmc
      use setout_gj_i, only: setout_gj
      use log_factorial_table, only: factt
      use mpi_utilities, only: mpix_shutdown
      use mpi_utilities, only: mpix_startup
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: input_status, nci, ncore_not_used, ncount1
      logical :: yes, use_args
      real(real64), pointer :: gjc_diag(:), dgjc_diag(:)
      character :: name*24
!-----------------------------------------------
!
      call mpix_startup(mpi_rank, mpi_process_count, processor_name, processor_name_length, ncount1, &
         "GJ_MPI", "This is the MPI Lande-factor program", &
         "isodata, descriptor Parquet/header, name.mix/name.cimix, name.wfn", &
         "name.gj or name.cgj", .true.)
!
      name = " "
      nci = 1
      ndef = 0
      use_args = .false.
!
      if (mpi_rank == 0) then
         use_args = command_argument_count() > 0
         if (use_args) then
            call property_state_arguments("gj_mpi", name, nci, input_status)
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
               write (stderr_unit, *) "GJ_MPI currently supports default settings only."
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
