

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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use memory_man
      use coun_c, only: count_context, count_trace_rank, double, thresh
      use def_c, only: nuclear_charge, speed_of_light_au, max_radial_points
      use default_c, only: double, inpci, ndef, ndump
      use grid_c
      use left_c, only: set
      use npar_c, only: double, nuclear_model_parameter_count, nuclear_model_parameters
      use orb_c
      use wave_c, only: radial_large_components, radial_small_components
      use whfrom_c, only: source
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i, only: getyn
      use tfpot_i, only: tfpot
      use prtrem_i, only: prtrem
      use getrsl_i, only: getrsl
      use frmrwf_i, only: frmrwf
      use frmtfp_i, only: frmtfp
      use frmhyd_i, only: frmhyd
      use summry_i, only: summry
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer  :: index(max_orbitals)
      integer :: j, i, k, nradial, nsubs, loc
      real(real64) :: con, fkk
      logical :: all, yes, existed
      logical, save :: modify = .false.
      character :: infile*128
!-----------------------------------------------

!
!   Set the threshold for node counting
!
      thresh = 0.05_real64
!
!   Set up the Thomas-Fermi potential
!
      call tfpot
!
!   Allocate storage to the arrays that store the subshell
!   radial wavefunction arrays; initialise these and all
!   associated arrays
!
      call alloc (radial_large_components,max_radial_points,orbital_count,"PF", "GENRWF")
      call alloc (radial_small_components,max_radial_points,orbital_count,"QF", "GENRWF")

      con = nuclear_charge/speed_of_light_au
      con = con*con

      do j = 1, orbital_count
         set(j) = .false.
         source(j) = " "
         radial_large_components(:radial_point_count,j) = 0.0_real64
         radial_small_components(:radial_point_count,j) = 0.0_real64

         k = abs(orbital_kappa(j))
         if (nuclear_model_parameter_count > 0) then
            orbital_origin_exponents(j) = dble(k)
         else if (nuclear_model_parameter_count == 0) then
            fkk = dble(k*k)
            if (fkk >= con) then
               orbital_origin_exponents(j) = sqrt(fkk - con)
            else
               write (stderr_unit, *) "LODRWF: Imaginary gamma parameter ", "for ", &
                  orbital_principal_n(j), orbital_symmetry_labels(j), " orbital;"
               write (stderr_unit, *) "the point model for the nucleus ", &
                  "is inappropriate for Z > ", speed_of_light_au, "."
               stop
            end if
         end if
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
         write (stderr_unit, *) "Read subshell radial wavefunctions. ", &
            "Choose one below"
         write (stderr_unit, *) "    1 -- GRASP92 File"
         write (stderr_unit, *) "    2 -- Thomas-Fermi"
         write (stderr_unit, *) "    3 -- Screened Hydrogenic"
         write (stderr_unit, *) "    4 -- Screened Hydrogenic [custom Z]"

         read (stdin_unit, *) nradial
         if (nradial<1 .or. nradial>4) then
            write (stderr_unit, *) nradial, "is not a valid choice, redo"
            goto 234
         end if

         if (nradial == 1) then
  345       continue
            write (stderr_unit, *) 'Enter the file name (Null then "rwfn.out")'
            read (stdin_unit, "(A)") infile
            if (len_trim(infile) == 0) infile = "rwfn.out"

            inquire(file=infile, exist=existed)
            if (.not.existed) then
               write (stderr_unit, *) ' File "', infile(1:len_trim(infile)), &
                  '" does not exist, redo'
               goto 345
            end if
         end if

         write (stderr_unit, *) "Enter the list of relativistic subshells:"
         open(unit=734, file="tmp_734", status="NEW")
         call getrsl (index, nsubs)
         close(734, status="DELETE")

         if (nradial == 1) then
            call frmrwf (index, nsubs, infile)
         else if (nradial == 2) then
            call frmtfp (index, nsubs)
         else if (nradial == 3) then
            call frmhyd (index, nsubs, modify,nradial)
         else if (nradial == 4) then
            call frmhyd (index, nsubs, modify,nradial)
         end if

         call prtrem (all)
         if (.not.all) goto 234

      end if
!
! All read. Let know, and allow modifying if non default
!
      write (stderr_unit, *) "All required subshell radial wavefunctions ", &
         " have been estimated:"
      call summry (stderr_unit)

      if (ndef == 0) then
         modify = .false.
      else
         write (stderr_unit, *) "Revise any of these estimates?"
         modify = getyn()
      end if

      if (modify) then
  456    continue
         write (stderr_unit, *) "Enter the list of subshells whose radial ", &
            "wavefunctions are to be revised:"
         open(unit=734, file="tmp_734", status="NEW")
         call getrsl (index, nsubs)
         close(734, status="DELETE")
         if (nsubs == 0) goto 456
         do j = 1, nsubs
            loc = index(j)
            set(loc) = .false.
            radial_large_components(:radial_point_count,loc) = 0.0_real64
            radial_small_components(:radial_point_count,loc) = 0.0_real64
            source(loc) = " "
         end do
         all = .false.
         goto 123
      end if

      if (ndef /= 0) call summry (24)

      return
      end subroutine genrwf
