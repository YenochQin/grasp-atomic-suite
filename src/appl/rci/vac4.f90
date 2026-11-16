!***********************************************************************
!                                                                      *
      subroutine vac4
!                                                                      *
!   This routine sets up the fourth-order vacuum polarization poten-   *
!   tial using equations (11) and (12) of L Wayne Fullerton and  G A   *
!   Rinker, Jr,  Phys  Rev  A 13 (1976) 1283-1287. The potential  is   *
!   accumulated in array  TC(I), I = 1, ..., N. It is transferred to   *
!   array TA in COMMON block TATB.                                     *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!               [RCI92]: FUNL.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 15 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use parameter_def,   only: max_radial_points
      use debug_c
      use def_c
      use grid_c, only: radial_grid, radial_point_count
      use npar_c
      use ncdist_c, only: zdist
      use tatb_c, only: mtp, ta, tb
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use funl_i
      use quad_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, nb2, nrows, ii, ii1, ii2
      real(double), dimension(max_radial_points) :: tc
      real(double) :: epsi, twocv, factor, ri, x, tci, rk, xk, xi, xm, xp
!-----------------------------------------------
!
!   Overall initialization
!
      epsi = precis*precis
      twocv = vacuum_speed_of_light_au + vacuum_speed_of_light_au
!
!   Potential for point nucleus: equation (12)
!
      factor = -nuclear_charge/(pi*vacuum_speed_of_light_au)**2
!
      tc(1) = 0.0d00
!
      i = 1
    1 continue
      i = i + 1
      ri = radial_grid(i)
      x = twocv*ri
      tci = (factor/ri)*funl(x,1)
      if (abs(tci) >= epsi) then
         tc(i) = tci
         if (i < radial_point_count) go to 1
      else
         tc(i:radial_point_count) = 0.0d00
      endif
!
!   Potential for finite nucleus: equation (11)
!
      if (nuclear_model_parameter_count == 2) then
!
         factor = -1.0d00/(pi*vacuum_speed_of_light_au**3)
!
         tc(1) = 0.0d00
!
         k = 1
    3    continue
         k = k + 1
!
         rk = radial_grid(k)
         xk = twocv*rk
         ta(1) = 0.0d00
!
         do i = 2, mtp
            xi = twocv*radial_grid(i)
            xm = abs(xk - xi)
            xp = xk + xi
            ta(i) = (funl(xm,0) - funl(xp,0))*zdist(i)
         end do
!
         call quad (x)
!
         x = x*factor/rk
!
!   Get out of the loop if the asymptotic region has been reached
!
         if (abs(x) >= epsi) then
            if (abs((tc(k)-x)/x) > 1.0d-03) then
               tc(k) = x
               if (k < radial_point_count) go to 3
            endif
         endif
!
      endif
!
      if (ldbpr(8)) then
         write (99, 300)
         nb2 = radial_point_count/2
         if (2*nb2 == radial_point_count) then
            nrows = nb2
         else
            nrows = nb2 + 1
         endif
         do ii = 1, nrows
            ii1 = ii
            ii2 = ii1 + nrows
            if (ii2 <= radial_point_count) then
               write (99, 301) radial_grid(ii1), tb(ii1), tc(ii1), radial_grid(ii2), tb(ii2), tc(&
                  ii2)
            else if (ii1 <= radial_point_count) then
               write (99, 301) radial_grid(ii1), tb(ii1), tc(ii1)
            endif
         end do
      endif
!
!   Generate total vacuum-polarization potential
!
      tb(:radial_point_count) = tc(:radial_point_count) + tb(:radial_point_count)
!
      return
!
  300 format(/,/,/,' ++++++++++ VAC4 ++++++++++'/,/,2(&
         ' -------- r -------- ----- VV2 (r) -----',' ----- VV4 (r) -----'))
  301 format(1p,6(1x,1d19.12))
      return
!
      end subroutine vac4
