!***********************************************************************
!                                                                      *
      subroutine getosd(isofile,name)
!                                                                      *
!   Interactively determines the data governing the transition prob-   *
!   lem.                                                               *
!                                                                      *
!   Written by Per Jonsson                Last revision: June 1996     *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_orbitals, max_radial_points
      use def_c
      use decide_c
      use default_c
      use foparm_c
      use grid_c, only: transformed_grid_step, asymptotic_grid_step, radial_point_count, radial_grid_scale
      use npar_c
      use osc_c, only: ltc
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
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character  :: name(2)*24
      character(len=*) :: isofile
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      logical :: yes
      character :: cunits*4, answer
!-----------------------------------------------
!
!   Open, check, load data from, and close the  .iso  file
!
!GG      CALL SETISO ('isodata')
      call setiso (isofile)
!
!   Determine the physical effects specifications
!
      if (ndef /= 0) then
         write (6, *) 'The physical speed of light in'
         write (6, *) ' atomic units is', vacuum_speed_of_light_au, ';'
         write (6, *) ' revise this value?'
         yes = getyn()
         if (yes) then
            write (6, *) 'Enter the revised value:'
            read (5, *) speed_of_light_au
         else
            speed_of_light_au = vacuum_speed_of_light_au
         endif
      else
         speed_of_light_au = vacuum_speed_of_light_au
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
      else
         ltc(5) = .true.
      endif
!
!      WRITE (6, *) 'Sort transitions by energy?'
!      YES = GETYN()
!      IF (YES) LTC(6) = .TRUE.
!
      if (ndef /= 0) then
         write (6, *) 'Einstein A and B coefficients are'
         write (6, *) ' printed in SI units; use Hartree'
         write (6, *) ' atomic units instead?'
         yes = getyn()
         if (yes) ltc(7) = .true.
      else
         ltc(7) = .false.
      endif
!
!   Determine the parameters controlling the radial grid
!
      if (nuclear_model_parameter_count == 0) then
         radial_grid_scale = exp((-65.0_real64/16.0_real64))/nuclear_charge
         transformed_grid_step = 0.5_real64**4
         radial_point_count = min(220,max_radial_points)
      else
!CFF     .. should be Z-dependent
         radial_grid_scale = 2.0e-6_real64/nuclear_charge
         transformed_grid_step = 5.0e-2_real64
         radial_point_count = max_radial_points
      endif
      asymptotic_grid_step = 0.0_real64
      if (ndef /= 0) then
         write (6, *) 'The default radial grid parameters'
         write (6, *) ' for this case are:'
         write (6, *) ' RNT = ', radial_grid_scale, ';'
         write (6, *) ' H = ', transformed_grid_step, ';'
         write (6, *) ' HP = ', asymptotic_grid_step, ';'
         write (6, *) ' N = ', radial_point_count, ';'
         write (6, *) ' revise these values?'
         yes = getyn()
         if (yes) then
            write (6, *) 'Enter RNT:'
            read (5, *) radial_grid_scale
            write (6, *) 'Enter H:'
            read (5, *) transformed_grid_step
            write (6, *) 'Enter HP:'
            read (5, *) asymptotic_grid_step
            write (6, *) 'Enter N:'
            read (5, *) radial_point_count
         endif
      endif
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      accy = transformed_grid_step**6
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
