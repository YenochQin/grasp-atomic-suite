!***********************************************************************
!                                                                      *
      subroutine getsmd(name)
!                                                                      *
!   Interactively determines the data governing the SMS problem.       *
!                                                                      *
!   Call(s) to: [LIB92]: NUCPOT, RADGRD, SETQIC.                       *
!               [RCI92]: SETISO, SETRWF.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: max_radial_points
      use decide_c
      use def_c
      use default_c, only: double, inpci, ndef, ndump
      use foparm_c, only: iccut
      use grid_c
      use npar_c, only: double, nuclear_model_parameter_count, nuclear_model_parameters
      use npot_c, only: double, first_grid_index_beyond_fermi_c, max_radial_points, negative_radius_times_potential
      use orb_c
      use wave_c
      use wfac_c, only: double, wfact
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i, only: getyn
      use setiso_i, only: setiso
      use setqic_i, only: setqic
      use radgrd_i, only: radgrd
      use nucpot_i, only: nucpot
      use setrwfa_i, only: setrwfa
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character, intent(in) :: name*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: yes
!-----------------------------------------------
!
!
!   Open, check, load data from, and close the  .iso  file
!
      call setiso ("isodata")
!
!   Determine the physical effects specifications
!

      if (ndef /= 0) then
         write (6, *) "The physical speed of light in"
         write (6, *) " atomic units is", vacuum_speed_of_light_au, ";"
         write (6, *) " revise this value?"
         yes = getyn()
         if (yes) then
            write (6, *) "Enter the revised value:"
            read (5, *) speed_of_light_au
         else
            speed_of_light_au = vacuum_speed_of_light_au
         end if
      else
         speed_of_light_au = vacuum_speed_of_light_au
      end if
!
      if (ndef /= 0) then
!
         write (6, *) "Treat contributions of some CSFs"
         write (6, *) " as first-order perturbations?"
         yes = getyn()
         if (yes) then
            lfordr = .true.
            write (6, *) "The contribution of CSFs"
            write (6, *) " 1 -- ICCUT will be treated"
            write (6, *) " variationally; the remainder"
            write (6, *) " perturbatively; enter ICCUT:"
            read (5, *) iccut
         else
            lfordr = .false.
            iccut = 0
         end if
      else
         lfordr = .false.
         iccut = 0
      end if
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
      end if
      asymptotic_grid_step = 0.0_real64
      if (ndef /= 0) then
         write (6, *) "The default radial grid parameters"
         write (6, *) " for this case are:"
         write (6, *) " RNT = ", radial_grid_scale, ";"
         write (6, *) " H = ", transformed_grid_step, ";"
         write (6, *) " HP = ", asymptotic_grid_step, ";"
         write (6, *) " N = ", radial_point_count, ";"
         write (6, *) " revise these values?"
         yes = getyn()
         if (yes) then
            write (6, *) "Enter RNT:"
            read (5, *) radial_grid_scale
            write (6, *) "Enter H:"
            read (5, *) transformed_grid_step
            write (6, *) "Enter HP:"
            read (5, *) asymptotic_grid_step
            write (6, *) "Enter N:"
            read (5, *) radial_point_count
         end if
      end if
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      accy = transformed_grid_step**6
!
!   Set up the coefficients for the numerical procedures
!
      call setqic
!
!   Generate the radial grid and all associated arrays
!
      call radgrd
!
!   Generate $- r \times V_nuc (r)$
!
      call nucpot
!
!   Load the radial wavefunctions
!
!      CALL SETRWFA(NAME)
      call setrwfa (trim(name)//".w")
!
      return
      end subroutine getsmd
