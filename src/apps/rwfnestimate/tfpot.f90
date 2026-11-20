

!***********************************************************************
!                                                                      *
      subroutine tfpot
!                                                                      *
!   Calculation of the universal Thomas-Fermi potential.               *
!                                                                      *
!   Call(s) to: [LIB92]: DRAW.                                         *
!                                                                      *
!                                         Last revision: 09 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use atomic_calculation_state, only: nuclear_charge, electron_count
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use nuclear_potential_state, only: negative_radius_times_potential
      use scf_radial_potentials, only: yp
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use radial_ascii_plot, only: draw
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: radial_index, nb3, row_count, ii, ii1, ii2, ii3
      real(real64) :: third, wa, wb, wc, wd, we, wf
!-----------------------------------------------
!
!
      third = 1.0_real64/3.0_real64
!
      wa = nuclear_charge - dble(electron_count - 1)
      wb = max(nuclear_charge - wa,0.0_real64)
      wb = wb**third/0.8853_real64
      do radial_index = 1, radial_point_count
!
!   Rational function approximation to the universal Thomas-Fermi
!   function
!
         wc = sqrt(wb*radial_grid(radial_index))
         wd = wc*(0.60112_real64*wc + 1.81061_real64) + 1.0_real64
         we = wc*(wc*(wc*(wc*(0.04793_real64*wc + 0.21465_real64) + 0.77112_real64) + &
            1.39515_real64) + 1.81061_real64) + 1.0_real64
         wf = wd/we
         yp(radial_index) = (negative_radius_times_potential(radial_index)-wa)*wf*wf + wa
      end do
!
!   Debug printout
!
      if (ldbpr(26)) then
         write (99, '(/,/,/," Thomas-Fermi potential"/,/,3(" --------- r --------- ------ -r*V(r) ------"))')
         nb3 = radial_point_count/3
         if (3*nb3 == radial_point_count) then
            row_count = nb3
         else
            row_count = nb3 + 1
         end if
         do ii = 1, row_count
            ii1 = ii
            ii2 = ii1 + row_count
            ii3 = ii2 + row_count
            if (ii3 <= radial_point_count) then
               write (99, '(1p,6(1x,1d21.14))') &
                  radial_grid(ii1), yp(ii1), radial_grid(ii2), yp(ii2), radial_grid(ii3), yp(ii3&
                  )
            else if (ii2 <= radial_point_count) then
               write (99, '(1p,6(1x,1d21.14))') radial_grid(ii1), yp(ii1), radial_grid(ii2), yp(ii2)
            else
               write (99, '(1p,6(1x,1d21.14))') radial_grid(ii1), yp(ii1)
            end if
         end do
         call draw (yp, 1.0_real64, yp, 1.0_real64, radial_point_count)
      end if
!
      return
!
      return
!
      end subroutine tfpot
