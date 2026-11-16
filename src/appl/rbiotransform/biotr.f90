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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:17:22   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
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
      use genmcp_i
      use starttime_i, only: starttime
      use stoptime_i, only: stoptime
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
      integer :: ncount1
      real(real64) :: work(lwork1)
      real(real64) :: cishl(nlmax*nlmax)
      real(real64)   :: s(max_orbitals, max_orbitals)
      real(real64) :: cfshl(nlmax*nlmax)
      logical :: yes
      character :: name(2)*24
      character(len=128) :: isofile
!-----------------------------------------------
!
      isofile = 'isodata'
!   Debug flags
!
      ntestg = 0
      ntestl = 0
      ntest = max0(ntestl,ntestg)
!
      call starttime (ncount1, 'BIOTRANSFORM')
      write(*,*)
      write(*,*) 'RBIOTRANSFORM'
      write(*,*) 'This program transforms the initial and final wave'
      write(*,*) 'functions so that standard tensor albegra can be'
      write(*,*) 'used in evaluation of the transition parameters'
      write(*,*) 'Input files:  isodata, name1.c, name1.w, name1.(c)m'
      write(*,*) '              name2.c, name2.w, name2.(c)m         '
      print*,'              name1.TB, name2.TB (optional angular files)'
      write(*,*) 'Output files: name1.bw, name1.(c)bm, '
      write(*,*) '              name2.bw, name2.c(bm)'
      write(*,*) '              name1.TB, name2.TB (angular files)'
      write(*,*)

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
      call fname (name)
!
!   Open, check, load data from and close the initial state CSL file.
!
      call setcslb (name(1), ncore1,1)
!
!   Transfer the data to the initial state COMMON
!
      call tcsl (1)
!
!   Open, check, load data from and close the final state CSL file.
!
      call setcslb (name(2), ncore2,2)
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
      call lodrwfi (name(1), ntestg)
!
!   Read the the radial orbitals for the initial state
!
      call lodrwff (name(2), ntestg)
!
!   Calculate the radial overlap matrices
!
      write (*, *)
      write (*, *) ' ******************************************'
      write (*, *) '  Overlap matrix before orbital rotations'
      write (*, *) ' *****************************************'
      write (*, *)

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
      write (*, *)
      write (*, *) ' ****************************************'
      write (*, *) '  Overlap matrix after orbital rotations'
      write (*, *) ' ****************************************'
      write (*, *)

      call brkt
!
!  Write the transformed radial functions to file
!
      call radfile (name)
!
!   Obtain one-electron coupling coefficients for the initial state.
!   The coefficients are dumped on files one kappa in turn and
!   thus the different kappa can be manipulated independently.
!   The interface with the transformation part is in the routine mcp
!
      call genmcp (name(1), 1, ntestg, inpci)
!***** added by Yu Zou, Feb.18,2000 ***********
      do k = 1, kamax
         close(unit=80 + k, status='DELETE')
      end do
!***** added by Yu Zou, Feb.18,2000 ***********
!
!   Obtain one-electron coupling coefficients for the final state.
!   The coefficients are dumped on files one kappa in turn and
!   thus the different kappa can be manipulated independently.
!   The interface with the transformation part is in the routine mcp
!
      call genmcp (name(2), 2, ntestg, inpci)

      do k = 1, kamax
         close(unit=80 + k, status='DELETE')
      end do
      call stoptime (ncount1, 'BIOTRANSFORM')

      stop
      end program biotr
