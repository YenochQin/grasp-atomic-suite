!***********************************************************************
!                                                                      *
      subroutine radpar
!                                                                      *
!   This subroutine sets the parameters controlling the radial grid    *
!                                                                      *
!   Last revision: June 1996                                           *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_radial_points
      use def_c,           only: speed_of_light_au, vacuum_speed_of_light_au, nuclear_charge, accy
      use default_c, only: double, inpci, ndef, ndump
      use grid_c, only: double, radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use npar_c, only: double, nuclear_model_parameter_count, nuclear_model_parameters
      use mpi_c, only: double, mpi_bcast, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_integer, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i, only: getyn
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
       end if !myid=0
       call mpi_bcast(speed_of_light_au,1,mpi_double_precision,0,mpi_comm_world,mpi_error_code)
      else
         speed_of_light_au = vacuum_speed_of_light_au
      end if

      return
      end subroutine radpar
