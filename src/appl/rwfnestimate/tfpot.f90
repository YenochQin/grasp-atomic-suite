

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
      use debug_c
      use def_c, only: nuclear_charge, electron_count
      use grid_c
      use npot_c, only: negative_radius_times_potential
      use orb_c
      use pote_c, only: yp
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use draw_i
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, nb3, nrows, ii, ii1, ii2, ii3
      real(real64) :: third, wa, wb, wc, wd, we, wf
!-----------------------------------------------
!
!
      third = 1.0_real64/3.0_real64
!
      wa = nuclear_charge - dble(electron_count - 1)
      wb = max(nuclear_charge - wa,0.0_real64)
      wb = wb**third/0.8853_real64
      do i = 1, radial_point_count
!
!   Rational function approximation to the universal Thomas-Fermi
!   function
!
         wc = sqrt(wb*radial_grid(i))
         wd = wc*(0.60112_real64*wc + 1.81061_real64) + 1.0_real64
         we = wc*(wc*(wc*(wc*(0.04793_real64*wc + 0.21465_real64) + 0.77112_real64) + &
            1.39515_real64) + 1.81061_real64) + 1.0_real64
         wf = wd/we
         yp(i) = (negative_radius_times_potential(i)-wa)*wf*wf + wa
      end do
!
!   Debug printout
!
      if (ldbpr(26)) then
         write (99, 300)
         nb3 = radial_point_count/3
         if (3*nb3 == radial_point_count) then
            nrows = nb3
         else
            nrows = nb3 + 1
         endif
         do ii = 1, nrows
            ii1 = ii
            ii2 = ii1 + nrows
            ii3 = ii2 + nrows
            if (ii3 <= radial_point_count) then
               write (99, 301) radial_grid(ii1), yp(ii1), radial_grid(ii2), yp(ii2), radial_grid(ii3), yp(ii3&
                  )
            else if (ii2 <= radial_point_count) then
               write (99, 301) radial_grid(ii1), yp(ii1), radial_grid(ii2), yp(ii2)
            else
               write (99, 301) radial_grid(ii1), yp(ii1)
            endif
         end do
         call draw (yp, 1.0_real64, yp, 1.0_real64, radial_point_count)
      endif
!
      return
!
  300 format(/,/,/,' Thomas-Fermi potential'/,/,3(&
         ' --------- r --------- ------ -r*V(r) ------'))
  301 format(1p,6(1x,1d21.14))
      return
!
      end subroutine tfpot
