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
      use default_c, only: double, inpci, ndef, ndump
      use debug_c, only: ldbpr, cutoff
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i, only: getyn
      use setmc_i, only: setmc
      use setcon_i, only: setcon
      use fname_i, only: fname
      use mrgcsl_i, only: mrgcsl

      use getosd_i, only: getosd
      use strsum_i, only: strsum
      use factt_i, only: factt
      use oscl_i, only: oscl
      use starttime_i, only: starttime
      use stoptime_i, only: stoptime
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ntest, ncount1, ilbl
      logical :: yes
      character :: name(2)*24
      character(len=128) :: isofile
!-----------------------------------------------
!
      isofile = "isodata"
      ntest = 1001
      call starttime (ncount1, "RTRANSITION")
      write(*,*)
      write(*,*) "RTRANSITION"
      write(*,*) "This program computes transition parameters from"
      write(*,*) "transformed wave functions"
      write(*,*) "Input files:  isodata, name1_descriptors.parquet, name1.bw, name1.(c)bm"
      write(*,*) "              name2_descriptors.parquet, name2.bw, name2.(c)bm         "
      write(*,*) "              optional, name1.lsj.lbl, name2.lsj.lbl"
      write(*,*) "              name1.name2.KT (optional angular files)"
      write(*,*) "Output files: name1.name2.(c)t                       "
      write(*,*) "              optional, name1.name2.(c)t.lsj         "
      write(*,*) "              name1.name2.KT (angular files)         "
      write(*,*) "Here K is parity and rank of transition: -1,+1 etc   "




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
!Rasa -- start
      ldbpr = .false.
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
      call fname (name)
!
!   Open, check, load data from, and close, the initial and final state
!   .csl  files. These files are then merged to one file.
!
      call mrgcsl (name)
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
!   Get the remaining information
!
      call getosd (isofile,name)
!
!   Open and append a summary of the inputs to the  .sum  file
!
      ilbl = 0
      call strsum (name, inpci,ilbl)
!
!   Set up the table of logarithms of factorials
!
      call factt
!
!   Proceed with the transition calculation
!
      call oscl (name)
!
!   Print completion message
!
      call stoptime (ncount1, "RTRANSITION")
!
      stop
      end program bioscl
