

!***********************************************************************
!                                                                      *
      subroutine genrwf
!                                                                      *
!   Controls the computation of the subshell radial wavefunctions.     *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, GETRSL, GETYN.                         *
!               [ERWF]: FRMHYD, FRMRWF, FRMTFP, PRTREM, SUMMRY,        *
!                       TFPOT.                                         *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
!   Updated by Xinghong He                Last revision: 11 Nov 1997   *
!    for menu-driven and less questions
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use memory_man
      use coun_c
      use def_c, only: nuclear_charge, speed_of_light_au, max_radial_points
      use default_c
      use grid_c
      use left_c, only: set
      use npar_c
      use orb_c
      use wave_c, only: radial_large_components, radial_small_components
      use whfrom_c, only: source
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use tfpot_i
      use prtrem_i
      use getrsl_i
      use frmrwf_i
      use frmtfp_i
      use frmhyd_i
      use summry_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(max_orbitals) :: index
      integer :: j, i, k, nradial, nsubs, loc
      real(double) :: con, fkk
      logical :: all, yes, modify=.false., existed
      character :: infile*128
!-----------------------------------------------

!
!   Set the threshold for node counting
!
      thresh = 0.05
!
!   Set up the Thomas-Fermi potential
!
      call tfpot
!
!   Allocate storage to the arrays that store the subshell
!   radial wavefunction arrays; initialise these and all
!   associated arrays
!
      call alloc (radial_large_components,max_radial_points,orbital_count,'PF', 'GENRWF')
      call alloc (radial_small_components,max_radial_points,orbital_count,'QF', 'GENRWF')

      con = nuclear_charge/speed_of_light_au
      con = con*con

      do j = 1, orbital_count
         set(j) = .false.
         source(j) = ' '
         radial_large_components(:radial_point_count,j) = 0.d0
         radial_small_components(:radial_point_count,j) = 0.d0

         k = abs(orbital_kappa(j))
         if (nuclear_model_parameter_count > 0) then
            orbital_origin_exponents(j) = dble(k)
         else if (nuclear_model_parameter_count == 0) then
            fkk = dble(k*k)
            if (fkk >= con) then
               orbital_origin_exponents(j) = sqrt(fkk - con)
            else
               write (stderr_unit, *) 'LODRWF: Imaginary gamma parameter ', 'for ', &
                  orbital_principal_n(j), orbital_symmetry_labels(j), ' orbital;'
               write (stderr_unit, *) 'the point model for the nucleus ', &
                  'is inappropriate for Z > ', speed_of_light_au, '.'
               stop
            endif
         endif
      end do
!
!   Write out the complete list of subshell radial wave functions
!
      call prtrem (all)
!
! Direct to read radial functions till finish
!
  123 continue
      if (.not.all) then
  234    continue
         write (stderr_unit, *)
         write (stderr_unit, *) 'Read subshell radial wavefunctions. ', &
            'Choose one below'
         write (stderr_unit, *) '    1 -- GRASP92 File'
         write (stderr_unit, *) '    2 -- Thomas-Fermi'
         write (stderr_unit, *) '    3 -- Screened Hydrogenic'
         write (stderr_unit, *) '    4 -- Screened Hydrogenic [custom Z]'

         read (stdin_unit, *) nradial
         if (nradial<1 .or. nradial>4) then
            write (stderr_unit, *) nradial, 'is not a valid choice, redo'
            go to 234
         endif

         if (nradial == 1) then
  345       continue
            write (stderr_unit, *) 'Enter the file name (Null then "rwfn.out")'
            read (stdin_unit, '(A)') infile
            if (len_trim(infile) == 0) infile = 'rwfn.out'

            inquire(file=infile, exist=existed)
            if (.not.existed) then
               write (stderr_unit, *) ' File "', infile(1:len_trim(infile)), &
                  '" does not exist, redo'
               go to 345
            endif
         endif

         write (stderr_unit, *) 'Enter the list of relativistic subshells:'
         open(unit=734, file='tmp_734', status='NEW')
         call getrsl (index, nsubs)
         close(734, status='DELETE')

         if (nradial == 1) then
            call frmrwf (index, nsubs, infile)
         else if (nradial == 2) then
            call frmtfp (index, nsubs)
         else if (nradial == 3) then
            call frmhyd (index, nsubs, modify,nradial)
         else if (nradial == 4) then
            call frmhyd (index, nsubs, modify,nradial)
         endif

         call prtrem (all)
         if (.not.all) go to 234

      endif
!
! All read. Let know, and allow modifying if non default
!
      write (stderr_unit, *) 'All required subshell radial wavefunctions ', &
         ' have been estimated:'
      call summry (stderr_unit)

      if (ndef == 0) then
         modify = .false.
      else
         write (stderr_unit, *) 'Revise any of these estimates?'
         modify = getyn()
      endif

      if (modify) then
  456    continue
         write (stderr_unit, *) 'Enter the list of subshells whose radial ', &
            'wavefunctions are to be revised:'
         open(unit=734, file='tmp_734', status='NEW')
         call getrsl (index, nsubs)
         close(734, status='DELETE')
         if (nsubs == 0) go to 456
         do j = 1, nsubs
            loc = index(j)
            set(loc) = .false.
            radial_large_components(:radial_point_count,loc) = 0.d0
            radial_small_components(:radial_point_count,loc) = 0.d0
            source(loc) = ' '
         end do
         all = .false.
         go to 123
      endif

      if (ndef /= 0) call summry (24)

      return
      end subroutine genrwf
