!***********************************************************************
!                                                                      *
      subroutine nucpot
!                                                                      *
!   Evaluate the nuclear potential for point and Fermi models.         *
!                                                                      *
!   Call(s) to: [LIB92] ES.                                            *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 05 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:02   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use debug_c
      use def_c,           only: z, pi
      use grid_c
      use npar_c
      use npot_c,          only: zz, nnuc
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use es_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, nb3, nrows, ii, ii1, ii2, ii3
      real(double) :: c, a, abc, tabc, abc2, thabc2, abc3, cba, pi2, hpiac2, &
         h3php, s2mcba, s3mcba, sabc3, dmsas, en, zbn, ri, rmc, rmcba, rbc, &
         s2rcba, s3rcba
      logical :: set
!-----------------------------------------------
!
!
!   Point nucleus
!
      if (nparm == 0) then
!
         zz(:n) = z
!
!   Fermi distribution
!
      else if (nparm == 2) then
!
         c = parm(1)
         a = parm(2)
         abc = a/c
         tabc = 2.0d00*abc
         abc2 = abc*abc
         thabc2 = 3.0d00*abc2
         abc3 = abc2*abc
         cba = c/a
         pi2 = pi*pi
         hpiac2 = 0.5d00*pi2*abc2
         h3php = 1.5d00 + hpiac2
         call es ((-cba), s2mcba, s3mcba)
         sabc3 = 6.0d00*abc3
         dmsas = -sabc3*s3mcba
         en = 1.0d00 + abc2*pi2 + dmsas
         zbn = z/en
!
         set = .false.
         do i = 1, n
            ri = r(i)
            rmc = ri - c
            rmcba = rmc/a
            rbc = ri/c
            if (rbc <= 1.0d00) then
               call es (rmcba, s2rcba, s3rcba)
               zz(i) = zbn*(dmsas + sabc3*s3rcba + rbc*(h3php - thabc2*s2rcba&
                   - 0.5d00*rbc*rbc))
            else
               if (.not.set) then
                  nnuc = i
                  set = .true.
               endif
               call es ((-rmcba), s2rcba, s3rcba)
               zz(i) = z*(1.0d00 + thabc2*(rbc*s2rcba + tabc*s3rcba)/en)
            endif
         end do
      endif
!
      if (ldbpr(2)) then
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
               write (99, 301) r(ii1), zz(ii1), r(ii2), zz(ii2), r(ii3), zz(ii3&
                  )
            else if (ii2 <= n) then
               write (99, 301) r(ii1), zz(ii1), r(ii2), zz(ii2)
            else
               write (99, 301) r(ii1), zz(ii1)
            endif
         end do
      endif
!
      return
!
  300 format(/,'From SUBROUTINE NUCPOT:'/,3(&
         ' -------- r -------- ----- -r*V(r) -----'))
  301 format(1p,6(1x,1d19.12))
      return
!
      end subroutine nucpot
