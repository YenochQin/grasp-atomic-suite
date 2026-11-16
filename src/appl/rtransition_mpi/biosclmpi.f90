!***********************************************************************
!                                                                      *
      program bioscl
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use default_c
      use debug_c, only: ldbpr, cutoff
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use setmc_i
      use setcon_i
      use fname_i
      use mrgcsl_i
      use setcslm_i
      use getosd_i
      use strsum_i
      use factt_i
      use oscl_i
      use mpi_c

      use suite_system, only: exit => terminate_process
      use mpiu_i, only: startmpi2
      use mpiu_i, only: stopmpi2
      use setiso_i, only: setiso
      use sys_chdir_i, only: sys_chdir
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ntest,ncount1,ilbl,lenperm,lentmp, iii
      logical :: yes
      character :: name(2)*24
      character :: fullname(2)*128
      character(len=3)  ::  idstring
      character(len=128) :: isofile,startdir,permdir,tmpdir
!-----------------------------------------------
!
!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================
      call startmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length, ncount1,           &
                      startdir, permdir, tmpdir, 'RTRANSITION_MPI')
      write (idstring, '(I3.3)') mpi_rank
      print*, tmpdir , ' = tmpdir'

      lenperm = len_trim (permdir)
      lentmp = len_trim (tmpdir)

!=======================================================================
      isofile = trim(startdir)//'/isodata'
      ntest = 1001

      if (mpi_rank .eq. 0) then
         write (6, *)
         write (6, *) ' Default settings?'
         yes = getyn()
         write (6, *)
         if (yes) then
            ndef = 0
            ndump = 1
         else
            ndef = 1
            write (6, *) ' Dump angular data to file?'
            yes = getyn()
            if (yes) then
               ndump = 1
            else
               ndump = 0
             endif
         endif
         write (6, *)
         write (6, *) ' Input from a CI calculation?'
         yes = getyn()
         write (6, *)
         if (yes) then
            inpci = 0
         else
            inpci = 1
         endif
      endif   !myid=0
      call mpi_bcast (ndef,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (ndump,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (inpci,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
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
      call setmc
!
!   Set up the physical constants
!
      call setcon
!
!   Obtain the names of the initial and final state files
!
      if (mpi_rank .eq. 0) call fname(name)
      call mpi_bcast (name,48,mpi_character,0,mpi_comm_world,mpi_error_code)
      fullname(1) = trim(startdir) // '/' // trim(name(1))
      fullname(2) = trim(startdir) // '/' // trim(name(2))
!
!   Open, check, load data from, and close, the initial and final state
!   .csl  files. These files are then merged to one file.
!
!cjb print FULLNAME1= FULLNAME2=
      if (mpi_rank .eq. 0) then
        print*,"FULLNAME1=",fullname(1)
        print*,"FULLNAME2=",fullname(2)
      end if

      call mrgcsl (fullname)
!
!   Open, check, load data from, and close, the merged .csl  file
!
      call setcslm
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
      call setiso(isofile)

!
!   Get the remaining information
!
      call getosd (fullname)
!
!   Open and append a summary of the inputs to the  .sum  file
!
      if (mpi_rank == 0) then
          iii = len_trim(startdir)
          call sys_chdir(trim(startdir),iii,mpi_error_code)
          if (mpi_error_code.ne.0) then
              print *, "error changing dir!"
              call exit(1)
          endif
          call strsum(name,inpci,ilbl)
          iii = len_trim(tmpdir)
          call sys_chdir(trim(tmpdir)//'/'//idstring,iii+4,mpi_error_code)
          if (mpi_error_code.ne.0) then
              print *, "Error changing dir!"
              call exit(1)
          endif
      endif  !myid = 0
!
!   Set up the table of logarithms of factorials
!
      call factt
!
!   Proceed with the transition calculation
!
      call oscl (name,fullname,tmpdir,startdir,idstring)
!
!   Print completion message
!
      if (mpi_rank .eq. 0) then
          print *
          print *, 'RTRANSITION_MPI: Execution complete.'
      endif
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================

      call stopmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length,                     &
                           ncount1, 'RTRANSITION_MPI')
!=======================================================================
!
      stop
      end program bioscl
