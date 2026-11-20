

!***********************************************************************
!                                                                      *
      subroutine getinf
!                                                                      *
!   Interactively determines data useful for generating estimates of   *
!   the subshell radial wavefunctions.                                 *
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, NUCPOT, RADGRD, SETQIC.                *
!               [RCI92]: SETISO.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 21 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use atomic_calculation_state, only: accy, max_radial_points, nuclear_charge, speed_of_light_au, vacuum_speed_of_light_au
      use calculation_input_defaults, only: inpci, ndef, ndump
      use radial_grid_state, only: radial_grid_scale, transformed_grid_step, asymptotic_grid_step, radial_point_count
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yes_no_input, only: getyn
      use isotope_data_input, only: setiso
      use radial_stencil_weights, only: setqic
      use radial_grid_setup, only: radgrd
      use nuclear_coulomb_potential, only: nucpot
      implicit none (type, external)
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: yes
!-----------------------------------------------
!
!
!
!   Open, check, load data from, and close the  .iso  file
!
      call setiso ("isodata")
!
!   Set defaults
!
      speed_of_light_au = vacuum_speed_of_light_au
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
!
      if (ndef /= 0) then
         write (stderr_unit, *) "Change the default speed of light ", &
            "or radial grid parameters?"
!
         yes = getyn()
         if (yes) then
!
!   Modify the speed of light
!
            write (stderr_unit, *) "The physical speed of light in ", &
               "atomic units is", vacuum_speed_of_light_au, ";"
            write (stderr_unit, *) "revise this value?"
            yes = getyn()
            if (yes) then
               write (stderr_unit, *) "Enter the revised value:"
               read (5, *) speed_of_light_au
            end if
!
!   Modify the parameters controlling the radial grid
!
            write (stderr_unit, *) "The default radial grid parameters ", &
               "for this case are:"
            write (stderr_unit, *) " RNT = ", radial_grid_scale, ";"
            write (stderr_unit, *) " H = ", transformed_grid_step, ";"
            write (stderr_unit, *) " HP = ", asymptotic_grid_step, ";"
            write (stderr_unit, *) " N = ", radial_point_count, ";"
            write (stderr_unit, *) "revise these values?"
            yes = getyn()
            if (yes) then
               write (stderr_unit, *) "Enter RNT:"
               read (5, *) radial_grid_scale
               write (stderr_unit, *) "Enter H:"
               read (5, *) transformed_grid_step
               write (stderr_unit, *) "Enter HP:"
               read (5, *) asymptotic_grid_step
               write (stderr_unit, *) "Enter N:"
               read (5, *) radial_point_count
            end if
!
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
!   Generate $- r \times V_ (r)$
!
      call nucpot
!
      return
      end subroutine getinf
