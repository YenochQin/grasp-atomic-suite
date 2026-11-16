!***********************************************************************
!                                                                      *
      subroutine vac2
!                                                                      *
!   This routine sets up the second-order vacuum polarization poten-   *
!   tial using  equations (1) and (4) of  L Wayne Fullerton and  G A   *
!   Rinker, Jr,  Phys Rev A  13 (1976) 1283-1287.  The  potential is   *
!   accumulated  in  array  TB(I), I = 1, ..., N  which is in COMMON   *
!   block /TATB/ .                                                     *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!               [RCI92]: FUNK.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use def_c
      use grid_c
      use npar_c
      use ncdist_c
      use tatb_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use funk_i
      use quad_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k
      real(double) :: epsi, twocv, factor, ri, x, tbi, rk, xk, xi, xm, xp
!-----------------------------------------------
!
!   Overall initialization
!
      epsi = precis*precis
      twocv = vacuum_speed_of_light_au + vacuum_speed_of_light_au
!
!   Potential for a point nucleus: equation (1)
!   (this is also the asymptotoc form for a finite nucleus)
!
      factor = -(2.0d00*nuclear_charge)/(3.0d00*pi*vacuum_speed_of_light_au)
!
      tb(1) = 0.0d00
!
      i = 1
    1 continue
      i = i + 1
!
      ri = radial_grid(i)
      x = twocv*ri
      tbi = (factor/ri)*funk(x,1)
!
      if (abs(tbi) >= epsi) then
         tb(i) = tbi
         if (i < radial_point_count) go to 1
      else
         tb(i:radial_point_count) = 0.0d00
      endif
!
!   Potential for a finite nucleus: equation (4)
!
      if (nuclear_model_parameter_count == 2) then
!
         factor = -2.0d00/(3.0d00*vacuum_speed_of_light_au**2)
!
!   Set up integrand
!
         tb(1) = 0.0d00
!
         k = 1
    3    continue
         k = k + 1
!
         rk = radial_grid(k)
         xk = twocv*rk
!
         ta(1) = 0.0d00
         do i = 2, mtp
            xi = twocv*radial_grid(i)
            xm = abs(xk - xi)
            xp = xk + xi
            ta(i) = (funk(xm,0) - funk(xp,0))*zdist(i)
         end do
!
         call quad (x)
!
         x = x*factor/rk
!
!   Get out of loop if the asymptotic value has been attained
!
         if (abs(x) >= epsi) then
            if (abs((x - tb(k))/x) > 1.0d-05) then
               tb(k) = x
               if (k < radial_point_count) go to 3
            endif
         endif
!
      endif
!
      return
      end subroutine vac2
