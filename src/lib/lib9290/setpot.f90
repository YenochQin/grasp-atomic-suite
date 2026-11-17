module setpot_i
    implicit none (type, external)
    private
    public :: setpot
contains
!***********************************************************************
!                                                                      *
      subroutine setpot(j, jp)
!                                                                      *
!   This  subroutine  sets  up the  arrays TF and TG for use by  the   *
!   subprograms IN, OUT, and SBSTEP.                                   *
!                                                                      *
!   The routine uses an algorithm that involves the orbital energy for *
!   finding a join.  When the energy estimate is such that the         *
!   criteria for the join is not met, the program sets the joining     *
!   point in the middle of the grid.                                   *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!      J:  (Input) Index of orbital                                    *
!      JP: (Output) Join point: point where TG changes sign            *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 16 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      use def_c, only: accy, double, max_radial_points, speed_of_light_au
      use grid_c, only: double, radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use int_c,           only: tf, tg
      use orb_c, only: double, orbital_energy_parameters
      use pote_c,          only: yp
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
      integer, intent(out) :: jp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(real64) :: dmhb2c, energy, enefac, rpi, yprpor
      logical :: notset
!-----------------------------------------------
!
!
      notset = .true.
!
!   Define constants
!
      dmhb2c = -transformed_grid_step/(2.0_real64*speed_of_light_au)
      energy = orbital_energy_parameters(j)
      enefac = 2.0_real64*speed_of_light_au*speed_of_light_au - energy
!
!   Set up arrays TF and TG
!
!   Since TF(1) and TG(1) are never used, set them
!   to some arbitrary value
!
      jp = 0
      tf(1) = 0.0_real64
      tg(1) = 0.0_real64
      do i = 2, radial_point_count
         rpi = radial_grid_derivative(i)
         yprpor = yp(i)*radial_derivative_over_radius(i)
         tf(i) = dmhb2c*(enefac*rpi + yprpor)
         tg(i) = dmhb2c*(energy*rpi - yprpor)
         if (.not.notset) cycle
         if (abs(sign(1.0_real64,tg(i))+sign(1.0_real64,tg(1))) >= accy) cycle
         jp = i
         notset = .false.
      end do
!
!   Trap for inappropriate grid
!
      if (jp == 0) then
         write (stderr_unit, *) "SETPOT: Join set to NNNP/2 = ", max_radial_points/2
         jp = max_radial_points/2
      end if
!
      return
      end subroutine setpot
end module setpot_i
