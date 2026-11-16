!***********************************************************************
!                                                                      *
      subroutine radgrd
!                                                                      *
!   This routine sets up the radial grid  R  and the associated arr-   *
!   ays  RP  and  RPOR  in the COMMON block  /GRID/. Different grids   *
!   are generated depending on whether or not  HP  is  zero .          *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:19   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use debug_c
      use def_c, only: precis
      use grid_c
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: np10, i, nb2, nrows, ii, ii1, ii2
      real(double) :: eph, ett, ettm1, epslon, a, rlast, rest, t, rests, fofr, &
         fpri, delr
!-----------------------------------------------
!
!
!   RPOR(1) is never used in the program: it is arbitrarily
!   set to zero
!
      np10 = n + 10
      r(1) = 0.0d00
      rpor(1) = 0.0d00
!
!   Now set up the grids
!
      if (hp == 0.0d00) then
!        default comes here
!
!   Exponential grid if HP is zero
!
!   Initializations
!
         rp(1) = rnt
         eph = exp(h)
         ett = 1.0d00
!
!   Set up the arrays R, RP, RPOR
!
         do i = 2, np10
            ett = eph*ett
            ettm1 = ett - 1.0d00
            r(i) = rnt*ettm1
            rp(i) = rnt*ett
            rpor(i) = ett/ettm1
         end do
!
      else
!
!   Asymptotically-linear exponential grid otherwise:
!
!   Initializations
!
         epslon = 1.0d03*precis
         a = h/hp
         rp(1) = rnt/(a*rnt + 1.0d00)
         rlast = 0.0d00
         rest = 0.0d00
!
!   Set up the arrays R, RP, RPOR
!
         do i = 2, np10
!
            t = h*dble(i - 1)
!
!   Solve the implicit equation for R using the Newton-Raphson
!   method
!
    2       continue
            rests = rest + rnt
            fofr = log(rests/rnt) + a*rest - t
            fpri = rests/(a*rests + 1.0d00)
            delr = -fofr*fpri
            rest = rlast + delr
!
            if (abs(delr/rest) < epslon) then
               r(i) = rest
               rests = rest + rnt
               fpri = rests/(a*rests + 1.0d00)
               rp(i) = fpri
               rpor(i) = fpri/rest
            else
               rlast = rest
               go to 2
            endif
!
         end do
!
      endif
!
!   Debug printout
!
      if (ldbpr(1)) then
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
               write (99, 301) r(ii1), rp(ii1), rpor(ii1), r(ii2), rp(ii2), &
                  rpor(ii2)
            else if (ii1 <= n) then
               write (99, 301) r(ii1), rp(ii1), rpor(ii1)
            endif
         end do
      endif
!
      return
!
  300 format(/,'From SUBROUTINE RADGRD:'/,2(&
         ' -------- r -------- -------- r'' -------',' ------- r''/r ------'))
  301 format(1p,6(1x,1d19.12))
      return
!
      end subroutine radgrd
