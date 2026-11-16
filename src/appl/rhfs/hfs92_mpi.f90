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
      use vast_kind_param, only: double
      use default_c
      use iounit_c
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use setdbg_i
      use setmc_i
      use setcon_i
      use setsum_i
      use setcsla_i
      use gethfd_i
      use getmixblock_i
      use strsum_i
      use factt_i
      use hfsgg_mpi_i
      use mpiu_i, only: gdmaxmpi_root
      use mpiu_i, only: mpix_shutdown
      use mpiu_i, only: mpix_startup
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: argc, iarg, k, nci, ncore_not_used, ncount1
      integer :: nstage1, nstage2, ncount_rate, ncount_max
      logical :: yes, use_args
      real(double) :: time_setcsla, time_gethfd, time_getmixblock, time_hfsgg
      real(double), dimension(1) :: timbuf
      character(len=24) :: arg
      character :: name*24
!-----------------------------------------------
!
      call mpix_startup(myid, nprocs, host, lenhost, ncount1, &
         'RHFS_MPI', 'This is the MPI hyperfine structure program', &
         'isodata, name.c, name.(c)m, name.w', &
         'name.(c)h, name.(c)hoffd', .true.)
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
                  write (istde, *) 'Usage: rhfs_mpi [name] [--ci|--nonci]'
                  write (istde, *) 'If no arguments are given, interactive mode is used.'
                  ndef = -1
               else if (iarg /= 1 .or. arg(1:2) == '--') then
                  write (istde, *) 'Unrecognized argument: ', trim(arg)
                  write (istde, *) 'Usage: rhfs_mpi [name] [--ci|--nonci]'
                  ndef = -1
               endif
            end do
            if (ndef == 0 .and. len_trim(name) == 0) then
               write (istde, *) 'Usage: rhfs_mpi [name] [--ci|--nonci]'
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
               write (istde, *) 'RHFS_MPI currently supports default settings only.'
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
      if (myid == 0) call setsum(name, nci)
      call system_clock(nstage1, ncount_rate, ncount_max)
      call setcsla(name, ncore_not_used)
      call system_clock(nstage2, ncount_rate, ncount_max)
      time_setcsla = dble(nstage2 - nstage1) / dble(ncount_rate)
      timbuf(1) = time_setcsla
      call gdmaxmpi_root(timbuf, 1)
      if (myid == 0) time_setcsla = timbuf(1)
!
      call system_clock(nstage1, ncount_rate, ncount_max)
      call gethfd(name)
      call system_clock(nstage2, ncount_rate, ncount_max)
      time_gethfd = dble(nstage2 - nstage1) / dble(ncount_rate)
      timbuf(1) = time_gethfd
      call gdmaxmpi_root(timbuf, 1)
      if (myid == 0) time_gethfd = timbuf(1)
!
      call system_clock(nstage1, ncount_rate, ncount_max)
      call getmixblock(name, nci)
      call system_clock(nstage2, ncount_rate, ncount_max)
      time_getmixblock = dble(nstage2 - nstage1) / dble(ncount_rate)
      timbuf(1) = time_getmixblock
      call gdmaxmpi_root(timbuf, 1)
      if (myid == 0) time_getmixblock = timbuf(1)
      if (myid == 0) call strsum
      call factt
!
      call system_clock(nstage1, ncount_rate, ncount_max)
      call hfsgg_mpi
      call system_clock(nstage2, ncount_rate, ncount_max)
      time_hfsgg = dble(nstage2 - nstage1) / dble(ncount_rate)
      timbuf(1) = time_hfsgg
      call gdmaxmpi_root(timbuf, 1)
      if (myid == 0) time_hfsgg = timbuf(1)
!
      if (myid == 0) then
         write (6, '(A, F10.3, A)') 'RHFS_MPI setcsla wall time  (max rank): ', &
            time_setcsla, ' s'
         write (6, '(A, F10.3, A)') 'RHFS_MPI gethfd wall time   (max rank): ', &
            time_gethfd, ' s'
         write (6, '(A, F10.3, A)') 'RHFS_MPI getmix wall time   (max rank): ', &
            time_getmixblock, ' s'
         write (6, '(A, F10.3, A)') 'RHFS_MPI hfsgg wall time    (max rank): ', &
            time_hfsgg, ' s'
      endif
!
      call mpix_shutdown(myid, ncount1, 'RHFS_MPI')
      stop
      end program hfs92_mpi
