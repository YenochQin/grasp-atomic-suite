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
      use vast_kind_param, only: double
      use default_c
      use iounit_c
      use memory_man
      use mpi_c
      use prnt_c,          only: nvec
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use compute_gj_only_mpi_i
      use gethfd_i
      use getmixblock_i
      use getyn_i
      use setcon_i
      use setcsla_i
      use setdbg_i
      use setmc_i
      use setout_gj_i
      use factt_i
      use mpiu_i, only: mpix_shutdown
      use mpiu_i, only: mpix_startup
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: argc, iarg, k, nci, ncore_not_used, ncount1
      logical :: yes, use_args
      real(double), dimension(:), pointer :: gjc_diag, dgjc_diag
      character(len=24) :: arg
      character :: name*24
!-----------------------------------------------
!
      call mpix_startup(myid, nprocs, host, lenhost, ncount1, &
         'GJ_MPI', 'This is the MPI Lande-factor program', &
         'isodata, name.c, name.(c)m, name.w', &
         'name.gj or name.cgj', .true.)
!
      name = ' '
      nci = 1
      ndef = 0
      use_args = .false.
!
      if (myid == 0) then
         argc = command_argument_count()
         if (argc > 0) then
            use_args = .true.
            do iarg = 1, argc
               arg = ' '
               call get_command_argument(iarg, arg)
               if (iarg == 1 .and. arg(1:2) /= '--') then
                  name = arg
               else if (trim(arg) == '--ci') then
                  nci = 0
               else if (trim(arg) == '--nonci') then
                  nci = 1
               else if (trim(arg) == '--help' .or. trim(arg) == '-h') then
                  write (istde, *) 'Usage: gj_mpi [name] [--ci|--nonci]'
                  write (istde, *) 'If no arguments are given, interactive mode is used.'
                  ndef = -1
               else if (iarg /= 1 .or. arg(1:2) == '--') then
                  write (istde, *) 'Unrecognized argument: ', trim(arg)
                  write (istde, *) 'Usage: gj_mpi [name] [--ci|--nonci]'
                  ndef = -1
               endif
            end do
            if (ndef == 0 .and. len_trim(name) == 0) then
               write (istde, *) 'Usage: gj_mpi [name] [--ci|--nonci]'
               ndef = -1
            endif
         endif
!
         if (.not.use_args .and. ndef == 0) then
            write (istde, *)
            write (istde, *) 'Default settings?'
            yes = getyn()
            write (istde, *)
            if (yes) then
               ndef = 0
            else
               ndef = 1
               write (istde, *) 'GJ_MPI currently supports default settings only.'
            endif
         endif
!
         if (ndef == 0 .and. .not.use_args) then
   10       continue
            write (istde, *) 'Name of state'
            read (*, '(A)') name
            k = index(name,' ')
            if (k == 1) then
               write (istde, *) 'Names may not start with a blank'
               go to 10
            endif
!
            write (istde, *)
            write (istde, *) 'Mixing coefficients from a CI calc.?'
            yes = getyn()
            if (yes) then
               nci = 0
            else
               nci = 1
            endif
         endif
      endif
!
      call mpi_bcast(ndef, 1, mpi_integer, 0, mpi_comm_world, ierr)
      if (ndef /= 0) then
         call mpi_finalize(ierr)
         stop
      endif
!
      call mpi_bcast(name, len(name), mpi_character, 0, mpi_comm_world, ierr)
      call mpi_bcast(nci, 1, mpi_integer, 0, mpi_comm_world, ierr)
!
      call setdbg
      call setmc
      call setcon
      call setcsla(name, ncore_not_used)
      call gethfd(name)
      call getmixblock(name, nci)
      call factt
!
      call alloc(gjc_diag, nvec, 'GJC_DIAG', 'GJ_MPI')
      call alloc(dgjc_diag, nvec, 'DGJC_DIAG', 'GJ_MPI')
!
      call compute_gj_only_mpi(gjc_diag, dgjc_diag)
      if (myid == 0) call setout_gj(name, nci, gjc_diag, dgjc_diag)
!
      call dalloc(gjc_diag, 'GJC_DIAG', 'GJ_MPI')
      call dalloc(dgjc_diag, 'DGJC_DIAG', 'GJ_MPI')
!
      call mpix_shutdown(myid, ncount1, 'GJ_MPI')
      stop
      end program gj_mpi
