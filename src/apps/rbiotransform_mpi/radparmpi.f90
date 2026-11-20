!***********************************************************************
!                                                                      *
      subroutine radpar
      use, intrinsic :: iso_fortran_env, only: input_unit, output_unit
!                                                                      *
!   This subroutine sets the parameters controlling the radial grid    *
!                                                                      *
!   Last revision: June 1996                                           *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_radial_points
      use atomic_calculation_state,           only: speed_of_light_au, vacuum_speed_of_light_au, nuclear_charge, accy
      use calculation_input_defaults, only: inpci, ndef, ndump
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_integer, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yes_no_input, only: getyn
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical   :: yes
      character :: answer
!-----------------------------------------------
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
        if (mpi_rank == 0) then
         write (output_unit, *) "The default radial grid parameters"
         write (output_unit, *) " for this case are:"
         write (output_unit, *) " RNT = ", radial_grid_scale, ";"
         write (output_unit, *) " H = ", transformed_grid_step, ";"
         write (output_unit, *) " HP = ", asymptotic_grid_step, ";"
         write (output_unit, *) " N = ", radial_point_count, ";"
         write (output_unit, *) " revise these values?"
         yes = getyn()
         if (yes) then
            write (output_unit, *) "Enter RNT:"
            read (input_unit, *) radial_grid_scale
            write (output_unit, *) "Enter H:"
            read (input_unit, *) transformed_grid_step
            write (output_unit, *) "Enter HP:"
            read (input_unit, *) asymptotic_grid_step
            write (output_unit, *) "Enter N:"
            read (input_unit, *) radial_point_count
         end if
        end if !myid=0
        call mpi_bcast(radial_grid_scale,1,mpi_double_precision,0,mpi_comm_world,mpi_error_code)
        call mpi_bcast(transformed_grid_step,1,mpi_double_precision,0,mpi_comm_world,mpi_error_code)
        call mpi_bcast(asymptotic_grid_step,1,mpi_double_precision,0,mpi_comm_world,mpi_error_code)
        call mpi_bcast (radial_point_count, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
      end if
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      accy = transformed_grid_step**6
!

      if (ndef /= 0) then
       if(mpi_rank == 0) then
         write (output_unit, *) "The physical speed of light in"
         write (output_unit, *) " atomic units is", vacuum_speed_of_light_au, ";"
         write (output_unit, *) " revise this value?"
         yes = getyn()
         if (yes) then
            write (output_unit, *) "Enter the revised value:"
            read (input_unit, *) speed_of_light_au
         else
            speed_of_light_au = vacuum_speed_of_light_au
         end if
       end if !myid=0
       call mpi_bcast(speed_of_light_au,1,mpi_double_precision,0,mpi_comm_world,mpi_error_code)
      else
         speed_of_light_au = vacuum_speed_of_light_au
      end if

      return
      end subroutine radpar
