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
      use parameter_def,   only: nnnp
      use debug_c
      use def_c
      use grid_c, only: r, n
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
      real(double), dimension(nnnp) :: tc
      real(double) :: epsi, twocv, factor, ri, x, tci, rk, xk, xi, xm, xp
!-----------------------------------------------
!
!   Overall initialization
!
      epsi = precis*precis
      twocv = cvac + cvac
!
!   Potential for point nucleus: equation (12)
!
      factor = -z/(pi*cvac)**2
!
      tc(1) = 0.0d00
!
      i = 1
    1 continue
      i = i + 1
      ri = r(i)
      x = twocv*ri
      tci = (factor/ri)*funl(x,1)
      if (dabs(tci) >= epsi) then
         tc(i) = tci
         if (i < n) go to 1
      else
         tc(i:n) = 0.0d00
      endif
!
!   Potential for finite nucleus: equation (11)
!
      if (nparm == 2) then
!
         factor = -1.0d00/(pi*cvac**3)
!
         tc(1) = 0.0d00
!
         k = 1
    3    continue
         k = k + 1
!
         rk = r(k)
         xk = twocv*rk
         ta(1) = 0.0d00
!
         do i = 2, mtp
            xi = twocv*r(i)
            xm = dabs(xk - xi)
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
         if (dabs(x) >= epsi) then
            if (dabs((tc(k)-x)/x) > 1.0d-03) then
               tc(k) = x
               if (k < n) go to 3
            endif
         endif
!
      endif
!
      if (ldbpr(8)) then
         write (99, 300)
         nb2 = n/2
         if (2*nb2 == n) then
            nrows = nb2
         else
            nrows = nb2 + 1
         endif
         do ii = 1, nrows
            ii1 = ii
            ii2 = ii1 + nrows
            if (ii2 <= n) then
               write (99, 301) r(ii1), tb(ii1), tc(ii1), r(ii2), tb(ii2), tc(&
                  ii2)
            else if (ii1 <= n) then
               write (99, 301) r(ii1), tb(ii1), tc(ii1)
            endif
         end do
      endif
!
!   Generate total vacuum-polarization potential
!
      tb(:n) = tc(:n) + tb(:n)
!
      return
!
  300 format(/,/,/,' ++++++++++ VAC4 ++++++++++'/,/,2(&
         ' -------- r -------- ----- VV2 (r) -----',' ----- VV4 (r) -----'))
  301 format(1p,6(1x,1d19.12))
      return
!
      end subroutine vac4
