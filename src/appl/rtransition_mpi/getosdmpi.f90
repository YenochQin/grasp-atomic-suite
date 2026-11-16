!***********************************************************************
!                                                                      *
      subroutine getosd(name)
!                                                                      *
!   Interactively determines the data governing the transition prob-   *
!   lem.                                                               *
!                                                                      *
!   Written by Per Jonsson                Last revision: June 1996     *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: nnnw, nnnp
      use def_c
      use decide_c
      use default_c
      use foparm_c
      use grid_c, only: h, hp, n, rnt
      use npar_c
      use osc_c, only: ltc
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use setiso_i
      use getrmp_i
      use setqic_i
      use radgrd_i
      use lodrwfi_i
      use lodrwff_i
      use brkt_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character (len = 128), dimension(2) :: name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      logical :: yes
      character :: cunits*4, answer
!-----------------------------------------------
!
!   Determine the physical effects specifications
!
      if (ndef /= 0) then
         if (myid .eq. 0) then
            print *, 'The physical speed of light in'
            print *, ' atomic units is',cvac,';'
            print *, ' revise this value?'
            yes = getyn ()
            if (yes) then
               print *, 'Enter the revised value:'
               read *,c
            else
               c = cvac
            endif
         endif !myid=0
         call mpi_bcast(c,1,mpi_double_precision,0,mpi_comm_world,ierr)
      else
         c = cvac
      endif
!
      lfordr = .false.
      iccut = 0
!
!   Determine the multipolarity and parity of the transitions
!
      call getrmp
!
!   Determine the units for the printout and other options
!
      ltc = .false.
!
      if (ndef /= 0) then
       if (myid .eq.0) then
         write (6, *) 'Which units are to be used to'
         write (6, *) ' express the transition energies?'
         write (6, *) '    A    : Angstrom:'
         write (6, *) '    eV   : electron volts;'
         write (6, *) '    Hart : Hartree atomic units;'
         write (6, *) '    Hz   : Hertz;'
         write (6, *) '    Kays : Kaysers [cm**(-1)];'
    2    continue
         read (*, '(A)') cunits
         if (cunits(1:1) == 'A') then
            ltc(1) = .true.
         else if (cunits(1:2) == 'eV') then
            ltc(2) = .true.
         else if (cunits(1:4) == 'Hart') then
            ltc(3) = .true.
         else if (cunits(1:2) == 'Hz') then
            ltc(4) = .true.
         else if (cunits(1:4) == 'Kays') then
            ltc(5) = .true.
         else
            write (6, *) 'GETOSD: Unable to interpret string;'
            write (6, *) ' reenter ...'
            go to 2
         endif
       endif !myid=0
       call mpi_bcast(ltc(1),5,mpi_logical,0,mpi_comm_world,ierr)
      else
         ltc(5) = .true.
      endif
!
!      WRITE (6, *) 'Sort transitions by energy?'
!      YES = GETYN()
!      IF (YES) LTC(6) = .TRUE.
!
      if (ndef /= 0) then
       if (myid .eq. 0) then
         write (6, *) 'Einstein A and B coefficients are'
         write (6, *) ' printed in SI units; use Hartree'
         write (6, *) ' atomic units instead?'
         yes = getyn()
         if (yes) ltc(7) = .true.
       endif !myid=0
       call mpi_bcast(ltc(7),1,mpi_logical,0,mpi_comm_world,ierr)
      else
         ltc(7) = .false.
      endif
!
!   Determine the parameters controlling the radial grid
!
      if (nparm == 0) then
         rnt = exp((-65.0d00/16.0d00))/z
         h = 0.5d00**4
         n = min(220,nnnp)
      else
!CFF     ... should be Z-dependent
         rnt = 2.0d-06/z
         h = 5.0d-02
         n = nnnp
      endif
      hp = 0.0d00
      if (ndef /= 0) then
        if (myid .eq. 0) then
         write (6, *) 'The default radial grid parameters'
         write (6, *) ' for this case are:'
         write (6, *) ' RNT = ', rnt, ';'
         write (6, *) ' H = ', h, ';'
         write (6, *) ' HP = ', hp, ';'
         write (6, *) ' N = ', n, ';'
         write (6, *) ' revise these values?'
         yes = getyn()
         if (yes) then
            write (6, *) 'Enter RNT:'
            read (5, *) rnt
            write (6, *) 'Enter H:'
            read (5, *) h
            write (6, *) 'Enter HP:'
            read (5, *) hp
            write (6, *) 'Enter N:'
            read (5, *) n
         endif
        endif !myid=0
        call mpi_bcast(rnt,1,mpi_double_precision,0,mpi_comm_world,ierr)
        call mpi_bcast(h,1,mpi_double_precision,0,mpi_comm_world,ierr)
        call mpi_bcast(hp,1,mpi_double_precision,0,mpi_comm_world,ierr)
        call mpi_bcast (n, 1, mpi_integer, 0, mpi_comm_world, ierr)
      endif
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      accy = h**6
!
!   Set up the coefficients for the numerical proceduRes
!
      call setqic
!
!   Generate the radial grid and all associated arrays
!
      call radgrd
!
!   Load the initial state radial wavefunctions.
!
      call lodrwfi (name(1))
!
!   Load the final state radial wavefunctions.
!
      call lodrwff (name(2))

!   Construct the radial overlap matrix.
!
      call brkt

      return
      end subroutine getosd
