

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
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use debug_c
      use def_c, only: z, nelec
      use grid_c
      use npot_c, only: zz
      use orb_c
      use pote_c, only: yp
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use draw_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, nb3, nrows, ii, ii1, ii2, ii3
      real(double) :: third, wa, wb, wc, wd, we, wf
!-----------------------------------------------
!
!
      third = 1.0d00/3.0d00
!
      wa = z - dble(nelec - 1)
      wb = max(z - wa,0.0d00)
      wb = wb**third/0.8853d00
      do i = 1, n
!
!   Rational function approximation to the universal Thomas-Fermi
!   function
!
         wc = sqrt(wb*r(i))
         wd = wc*(0.60112d0*wc + 1.81061d0) + 1.0d00
         we = wc*(wc*(wc*(wc*(0.04793d0*wc + 0.21465d0) + 0.77112d0) + &
            1.39515d0) + 1.81061d0) + 1.0d00
         wf = wd/we
         yp(i) = (zz(i)-wa)*wf*wf + wa
      end do
!
!   Debug printout
!
      if (ldbpr(26)) then
         write (99, 300)
         nb3 = n/3
         if (3*nb3 == n) then
            nrows = nb3
         else
            nrows = nb3 + 1
         endif
         do ii = 1, nrows
            ii1 = ii
            ii2 = ii1 + nrows
            ii3 = ii2 + nrows
            if (ii3 <= n) then
               write (99, 301) r(ii1), yp(ii1), r(ii2), yp(ii2), r(ii3), yp(ii3&
                  )
            else if (ii2 <= n) then
               write (99, 301) r(ii1), yp(ii1), r(ii2), yp(ii2)
            else
               write (99, 301) r(ii1), yp(ii1)
            endif
         end do
         call draw (yp, 1.0d00, yp, 1.0d00, n)
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
