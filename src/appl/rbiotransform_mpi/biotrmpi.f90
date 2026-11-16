!***********************************************************************
!                                                                      *
!     This program transforms the initial and final state radial       *
!     orbitals to a biorthonormal set. The CI coefficients are         *
!     then counter transformed as to leave both the initial and        *
!     final state wave functions invariant                             *
!                                                                      *
!     General references:                                              *
!                                                                      *
!     P.A. Malmqvist, Int.J. of Quantum Chemistry, XXX, 479-94 (1986)  *
!     J. Olsen, M.R. Godefroid, P. Jonsson, P.A. Malmqvist and         *
!     C. Froese Fischer, Phys. Rev. E, 4499 (1995)                     *
!                                                                      *
      program biotr
!                                                                      *
!     Program written by                                               *
!                                                                      *
!     Per Jonsson, Department of Computer Science                      *
!                  Vanderbilt University, USA                          *
!                                                                      *
!     Date June 1996                                                   *
!                                                                      *
!     Modified by Gediminas Gaigalas for new spin-angular integration  *
!     and for reducing usage of CPU memory.        NIST, October 2017  *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_orbitals, max_radial_points
      use default_c,       only: ndef, ndump
      use sbdat_c,         only: nlmax, kamax, nshlii, nshlff
      use biorb_c,         only: nwii, nwff
      use wave_c,          only: pfii, qfii, pfff, qfff
      use cimat_c,         only: cici, cfci
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use setmc_i
      use setcon_i
      use setiso_i
      use radpar_i
      use radgrd_i
      use setqic_i
      use fname_i
      use setcslb_i
      use tcsl_i
      use kapdata_i
      use lodrwfi_i
      use lodrwff_i
      use brkt_i
      use gets_i
      use biotr1_i
      use radfile_i
      use factt_i, only: factt
      use mcpin_i, only: mcpin
      use mpiu_i, only: startmpi2
      use mpiu_i, only: stopmpi2
      use setcsla_i, only: setcsla
      use sys_chdir_i, only: sys_chdir
      use mcpout_i, only: mcpout
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!GG      INTEGER, PARAMETER :: LWORK1 = 100000
      integer, parameter :: lwork1 = 10000000
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ninshli(nlmax), ninshlf(nlmax)
      integer :: ntestg, ntestl, ntest, inpci, ncore1, ncore2, mxl, k
      integer :: ncount1,lenperm,lentmp,iii,ncore_not_used
      real(real64) :: work(lwork1)
      real(real64) :: cishl(nlmax*nlmax)
      real(real64)   :: s(max_orbitals, max_orbitals)
      real(real64) :: cfshl(nlmax*nlmax)
      logical :: yes
      character :: name(2)*24
      character(len=3)  ::  idstring
      character(len=128) :: isofile
      character(len=128) :: fullname(2)
      character (len = 128) :: startdir,permdir,tmpdir
!-----------------------------------------------
!
!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================
      call startmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length, ncount1,            &
                        startdir, permdir, tmpdir, 'RBIOTRANSFORM_MPI')
      write (idstring, '(I3.3)') mpi_rank
      print*, tmpdir , ' = tmpdir'

      isofile = trim(startdir)//'/isodata'
      lenperm = len_trim (permdir)
      lentmp = len_trim (tmpdir)

!=======================================================================

!   Debug flags
!
      ntestg = 0
      ntestl = 0
      ntest = max0(ntestl,ntestg)
!
      if (mpi_rank .eq. 0) then
      print *
      print *

      write (6, *) 'Default settings?'
      yes = getyn()
      write (6, *)
      if (yes) then
         ndef = 0
         ndump = 1
      else
         ndef = 1
         write (6, *) 'Dump angular data on file?'
         yes = getyn()
         write (6, *)
         if (yes) then
            ndump = 1
         else
            ndump = 0
         endif
      endif
      write (6, *) 'Input from a CI calculation?'
      yes = getyn()
      write (6, *)
      if (yes) then
         inpci = 0
      else
         inpci = 1
      endif
      endif !myid=0
      call mpi_bcast (ndef,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (ndump,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (inpci,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
!
!   Perform machine- and installation-dependent setup
!
      call setmc
!
!   Set up the physical constants
!
      call setcon
!
!   Open, check, load data from, and close the  .iso  file
!
      call setiso (isofile)
!
!   Determine the parameters controlling the radial grid
!
      call radpar
!
!   Generate the radial grid
!
      call radgrd
!
!   Set up the coefficients for the numerical procedures
!
      call setqic
!
!   Obtain the names of the initial and final state files
!   and open files where the transformed orbitals and CI
!   coefficients are to be dumped
!
      if (mpi_rank .eq. 0) call fname (name)
      call mpi_bcast (name,48,mpi_character,0,mpi_comm_world,mpi_error_code)
      print *, 'DBG: myid=',mpi_rank, ' name(1)=', name(1)
      print *, 'DBG: myid=',mpi_rank, ' name(2)=', name(2)
      fullname(1) = trim(startdir) // '/' // trim(name(1))
      fullname(2) = trim(startdir) // '/' // trim(name(2))
      print *, 'DBG: myid=',mpi_rank, ' startdir =', startdir
      print *, 'DBG: myid=',mpi_rank, ' fullname(1)=', fullname(1)
      print *, 'DBG: myid=',mpi_rank, ' fullname(2)=', fullname(2)
!
!   Open, check, load data from and close the initial state CSL file.
!
      call setcslb (fullname(1), ncore1,1)
!
!   Transfer the data to the initial state COMMON
!
      call tcsl (1)
!
!   Open, check, load data from and close the final state CSL file.
!
      call setcslb (fullname(2), ncore2,2)
!
!   Transfer the data to the final state COMMON
!
      call tcsl (2)
!
!   Determine the number of kappa quantum numbers and
!   the number of orbitals for each kappa quantum number
!   for the initial state and final states
!
      call kapdata (ntestg, ncore1, ncore2)
!
!   Read the the radial orbitals for the initial state
!
      call lodrwfi (fullname(1), ntestg)
!
!   Read the the radial orbitals for the initial state
!
      call lodrwff (fullname(2), ntestg)
!
!   Calculate the radial overlap matrices
!
      if (mpi_rank .eq. 0) then
      write (*, *)
      write (*, *) ' ******************************************'
      write (*, *) '  Overlap matrix before orbital rotations'
      write (*, *) ' *****************************************'
      write (*, *)
      endif!myid=0

      call brkt

      call gets (s, nwii, nwff)

!
!   Once we have the overlap matrices
!   we can manipulate the initial and final state separately.
!
      mxl = kamax
!
!. Calculate biorthonormal orbitals, and orbital matrix
!. for counter transformation of CI coefficients.
!
      call biotr1 (pfii, qfii, nshlii, ninshli, pfff, qfff, &
         nshlff, ninshlf, max_radial_points, kamax, work, lwork1, ntestg, &
         cishl, cici, cfshl, cfci)
      if (mpi_rank .eq. 0) then
      write (*, *)
      write (*, *) ' ****************************************'
      write (*, *) '  Overlap matrix after orbital rotations'
      write (*, *) ' ****************************************'
      write (*, *)
      endif !myid=0

      call brkt
!
!  Write the transformed radial functions to file
!
      call radfile(fullname)
!
!   Obtain one-electron coupling coefficients for the initial state.
!   The coefficients are dumped on files one kappa in turn and
!   thus the different kappa can be manipulated independently.
!   The interface with the transformation part is in the routine mcp
!
      iii = len_trim(startdir)
      call sys_chdir(trim(startdir),iii,mpi_error_code)
!GG      CALL SETCSLA(NAME(1), ncore_not_used)
      call setcsla(name(1), ncore_not_used,3)
      iii = len_trim(tmpdir)
      call sys_chdir(trim(tmpdir)//'/'//idstring,iii+4,mpi_error_code)
      call factt
      call mcpout(name(1),startdir,1,ntestg,inpci)
      call mcpin(name(1),startdir,1,ntestg,inpci)

!***** added by Yu Zou, Feb.18,2000 ***********
!R      DO K = 1,KAMAX
!R        close(UNIT=80+K,STATUS="DELETE")
!R      ENDDO
!***** added by Yu Zou, Feb.18,2000 ***********
!
!   Obtain one-electron coupling coefficients for the final state.
!   The coefficients are dumped on files one kappa in turn and
!   thus the different kappa can be manipulated independently.
!   The interface with the transformation part is in the routine mcp
!
      iii = len_trim(startdir)
      call sys_chdir(trim(startdir),iii,mpi_error_code)
!GG      CALL SETCSLA(NAME(2), ncore_not_used)
      call setcsla(name(2), ncore_not_used,3)
      iii = len_trim(tmpdir)
      call sys_chdir(trim(tmpdir)//'/'//idstring,iii+4,mpi_error_code)
      call mcpout(name(2),startdir,2,ntestg,inpci)
      call mcpin(name(2),startdir,2,ntestg,inpci)

!R      DO K = 1,KAMAX
!R        close(UNIT=80+K,STATUS="DELETE")
!R      ENDDO

!=======================================================================
!  Execution finished; Statistics output
!=======================================================================

      call stopmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length,             &
                           ncount1, 'RBIOTRANSFORM_MPI')
!=======================================================================
      stop
      end program biotr
