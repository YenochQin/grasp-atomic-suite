!***********************************************************************
!                                                                      *
      subroutine prwf(j)
!                                                                      *
!   Makes a (debug) printout of wave functions. There are two modes:   *
!                                                                      *
!      (1) J>0  - Used as a debug option in SOLVE, wavefunctions for   *
!                 orbital J are printed                                *
!      (2) J=0  - A  printout  of the grid and all wave functions is   *
!                 made                                                 *
!                                                                      *
!   Call(s) to: [LIB92]: DRAW.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: nnnp
      use def_c
      use grid_c
      use int_c
      use orb_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE draw_I
      use draw_i, only: draw
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nb2, nrows, ii, ii1, ii2, k, mfk, jgg
      real(double) :: cbz
!-----------------------------------------------
!
!
      cbz = c/z
!
      if (j > 0) then
!
!   Mode (1)
!
         write (99, 300) np(j), nh(j)
         nb2 = mtp0/2
         if (2*nb2 == mtp0) then
            nrows = nb2
         else
            nrows = nb2 + 1
         endif
         do ii = 1, nrows
            ii1 = ii
            ii2 = ii1 + nrows
            if (ii2 <= mtp0) then
               write (99, 301) r(ii1), p(ii1), q(ii1), r(ii2), p(ii2), q(ii2)
            else if (ii1 <= mtp0) then
               write (99, 301) r(ii1), p(ii1), q(ii1)
            endif
         end do
         call draw (p, 1.0d00, q, cbz, mtp0)
!
      else
!
!   Mode (2)
!
         do k = 1, nw
            write (99, 300) np(k), nh(k)
            mfk = mf(k)
            nb2 = mfk/2
            if (2*nb2 == mfk) then
               nrows = nb2
            else
               nrows = nb2 + 1
            endif
            do ii = 1, nrows
               ii1 = ii
               ii2 = ii1 + nrows
               if (ii2 <= mfk) then
                  write (99, 301) r(ii1), pf(ii1,k), qf(ii1,k), r(ii2), pf(ii2,&
                     k), qf(ii2,k)
               else if (ii1 <= mfk) then
                  write (99, 301) r(ii1), pf(ii1,k), qf(ii1,k)
               endif
            end do
            call draw (pf(:nnnp,k), 1.0d00, qf(:nnnp,k), cbz, mf(k))
!
         end do
!
      endif
!
      return
!
  300 format('1',1i2,1a2,' orbital:'/,2(&
         ' --------- r --------- ------- P (r) -------',&
         ' ------- Q (r) -------'))
  301 format(1p,6(1x,1d21.14))
      return
!
      end subroutine prwf
