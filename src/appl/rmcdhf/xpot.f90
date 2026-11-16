!***********************************************************************
!                                                                      *
      subroutine xpot(j)
!                                                                      *
!   This subroutine tabulates the exchange terms (the first terms on   *
!   the  right-hand  sides of eqs (14), I P Grant, B J McKenzie, P H   *
!   Norrington, D F Mayers, and N C Pyper, Computer  Phys  Commun 21   *
!   (1980) 211) for orbital J. The exchange terms are stored  in the   *
!   common arrays XP and XQ.                                           *
!                                                                      *
!   Call(s) to: [LIB92]: DRAW, YZK.                                    *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 10 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:59:40   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: keyorb
      use debug_c
      use def_c
      use grid_c
      use orb_c
      use pote_c
      use scf_c
      use tatb_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yzk_i
      use draw_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: j
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = keyorb
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, index, label, k, ioy1, ioy2, iorb, nb2, nrows, ii, ii1, ii2
      real(double) :: coeff, ctb
!-----------------------------------------------
!
!   Debug printout: header
!
      if (ldbpr(27) .or. ldbpr(28)) write (99, 300) np(j), nh(j)
!
!   Clear for accumulation of sums
!
      xp(:n) = 0.d0
      xq(:n) = 0.d0
!
!   Add contributions from exchange terms
!
      do index = 1, nxcof

         ! Decode information in label
         label = nxa(index)
         k = mod(label,key)
         label = label/key
         ioy1 = mod(label,key)
         label = label/key
         ioy2 = mod(label,key)
         iorb = label/key
         coeff = xa(index)

         ! Debug printout: composition
         if (ldbpr(27)) write (99, 301) k, coeff, np(ioy1), nh(ioy1), np(ioy2)&
            , nh(ioy2), np(iorb), nh(iorb)

         call yzk (k, ioy1, ioy2)
!
!   Accumulate contributions
!
         coeff = coeff/c
         !DO I = 1, MF(IORB)
         do i = 1, n
            ctb = coeff*tb(i)
            xp(i) = xp(i) + ctb*qf(i,iorb)
            xq(i) = xq(i) - ctb*pf(i,iorb)
         end do
      end do
!
!   Debug printout: potential functions
!
      if (ldbpr(28)) then
         write (99, 302)
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
               write (99, 303) r(ii1), xp(ii1), xq(ii1), r(ii2), xp(ii2), xq(&
                  ii2)
            else if (ii1 <= n) then
               write (99, 303) r(ii1), xp(ii1), xq(ii1)
            endif
         end do
         call draw (xp, 1.0d00, xq, c, n)
      endif
!
      return
!
  300 format(/,/,' Exchange potential contributions (coefficients will ',&
         ' be divided by C) for ',1i2,1a2,' orbital :'/,/)
  301 format(/,25x,'(',1i2,')'/,1x,1p,d21.14,'* Y    (',1i2,1a2,',',1i2,1a2,&
         ') ','* P (',1i2,1a2,')')
  302 format(/,/,31x,'(P)',19x,'(Q)',41x,'(P)',19x,'(Q)'/,2(&
         ' --------- r --------- ------ X  (r) -------',&
         ' ------ X  (r) -------'))
  303 format(1p,6(1x,1d21.14))
      return
!
      end subroutine xpot
