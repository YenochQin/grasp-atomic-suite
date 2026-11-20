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
      use calculation_input_defaults, only: inpci, ndef, ndump
      use debug_output_settings, only: ldbpr, cutoff
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yes_no_input, only: getyn
      use machine_constants, only: setmc
      use physical_constants, only: setcon
      use fname_i, only: fname
      use mrgcsl_i, only: mrgcsl

      use getosd_i, only: getosd
      use strsum_i, only: strsum
      use log_factorial_table, only: factt
      use oscl_i, only: oscl
      use mpi_state, only: mpi_bcast, mpi_character, mpi_comm_world, mpi_error_code, mpi_integer, &
         & mpi_process_count, mpi_rank, processor_name, processor_name_length

      use posix_services, only: exit => terminate_process
      use mpi_utilities, only: startmpi2
      use mpi_utilities, only: stopmpi2
      use isotope_data_input, only: setiso
      use mpi_change_directory, only: sys_chdir
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ntest,start_clock_count,ilbl,lenperm,lentmp, iii
      logical :: yes
      character :: name(2)*24
      character :: fullname(2)*128
      character(len=3)  ::  block_label_text
      character(len=128) :: isotope_file_path,startup_directory,permdir,tmpdir
!-----------------------------------------------
!
!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================
      call startmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length, start_clock_count,           &
                      startup_directory, permdir, tmpdir, "RTRANSITION_MPI")
      write (block_label_text, "(I3.3)") mpi_rank
      print*, tmpdir , " = tmpdir"

      lenperm = len_trim (permdir)
      lentmp = len_trim (tmpdir)

!=======================================================================
      isotope_file_path = trim(startup_directory)//"/isodata"
      ntest = 1001

      if (mpi_rank == 0) then
         write (6, *)
         write (6, *) " Default settings?"
         yes = getyn()
         write (6, *)
         if (yes) then
            ndef = 0
            ndump = 1
         else
            ndef = 1
            write (6, *) " Dump angular data to file?"
            yes = getyn()
            if (yes) then
               ndump = 1
            else
               ndump = 0
             end if
         end if
         write (6, *)
         write (6, *) " Input from a CI calculation?"
         yes = getyn()
         write (6, *)
         if (yes) then
            inpci = 0
         else
            inpci = 1
         end if
      end if   !myid=0
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
      if (mpi_rank == 0) call fname(name)
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

      call mrgcsl (fullname)
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
      call setiso(isotope_file_path)

!
!   Get the remaining information
!
      call getosd (fullname)
!
!   Open and append a summary of the inputs to the  .sum  file
!
      if (mpi_rank == 0) then
          iii = len_trim(startup_directory)
          call sys_chdir(trim(startup_directory),iii,mpi_error_code)
          if (mpi_error_code/=0) then
              print *, "error changing dir!"
              call exit(1)
          end if
          ilbl = 0
          call strsum(name,inpci,ilbl)
          iii = len_trim(tmpdir)
          call sys_chdir(trim(tmpdir)//"/"//block_label_text,iii+4,mpi_error_code)
          if (mpi_error_code/=0) then
              print *, "Error changing dir!"
              call exit(1)
          end if
      end if  !myid = 0
!
!   Set up the table of logarithms of factorials
!
      call factt
!
!   Proceed with the transition calculation
!
      call oscl (name,fullname,tmpdir,startup_directory,block_label_text)
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

      call stopmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length,                     &
                           start_clock_count, "RTRANSITION_MPI")
!=======================================================================
!
      stop
      end program bioscl
