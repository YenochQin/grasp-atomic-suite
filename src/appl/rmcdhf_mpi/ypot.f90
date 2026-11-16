!***********************************************************************
!                                                                      *
      subroutine ypot(j)
!                                                                      *
!   This subroutine  tabulates the  potential function Y(r) (Eq (14)   *
!   in  I P Grant, B J McKenzie, P H Norrington, D F Mayers, and N C   *
!   Pyper, Computer  Phys  Commun  21 (1980) 211) for orbital J. The   *
!   function is tabulated in the COMMON array YP.                      *
!                                                                      *
!   Call(s) to: [LIB92]: DRAW, YZK                                     *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 05 Aug 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: keyorb
      use debug_c
      use grid_c
      use npot_c
      use orb_c
      use pote_c
      use scf_c
      use tatb_c
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yzk_i
      use draw_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = keyorb
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, index, label, k, ioy1, ioy2, nb3, nrows, ii, ii1, ii2, ii3
      real(double) :: coeff
!-----------------------------------------------
!   Debug printout: composition

      if (ldbpr(29) .or. ldbpr(30)) write (99, 300) np(j), nh(j)
!
!   Initialize array YP with the nuclear potential piece
!
!   Since YA() below  contains contributions from THIS node only,
!   the initialization should be in consistence with that.

      yp(:n) = zz(:n) / nprocs

      do index = 1, nycof

         ! Decode information in label
         label = nya(index)
         k = mod(label,key)
         label = label/key
         ioy1 = mod(label,key)
         ioy2 = label/key
         coeff = ya(index)

         if (ldbpr(29)) write (99, 301) k, coeff, np(ioy1), nh(ioy1), np(ioy2)&
            , nh(ioy2)

         call yzk (k, ioy1, ioy2)                ! Accumulate contributions
         yp(:n) = yp(:n) - coeff*tb(:n)
      end do
!
!   Debug printout
!
      if (ldbpr(30)) then
         write (99, 302)
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
               write (99, 303) r(ii1),yp(ii1),r(ii2),yp(ii2),r(ii3),yp(ii3)
            else if (ii2 <= n) then
               write (99, 303) r(ii1),yp(ii1),r(ii2),yp(ii2)
            else
               write (99, 303) r(ii1),yp(ii1)
            endif
         end do
         call draw (yp, 1.0d00, yp, 0.0d00, n)
      endif
!
      return
!
  300 format(/,/,' Direct potential for ',1i2,1a2,' orbital :'/,/)
  301 format(/,25x,'(',1i2,')'/,1x,1p,d21.14,'* Y    (',1i2,1a2,',',1i2,1a2,')'&
         )
  302 format(/,/,3(' --------- r --------- ------- Y (r) -------'))
  303 format(1p,6(1x,1d21.14))
      return
!
      end subroutine ypot
