module interp_i
    implicit none (type, external)
    private
    public :: interp
contains
!***********************************************************************
!                                                                      *
      subroutine interp(xarr, yarr, narr, xval, yval, accy)
!                                                                      *
!   This routine returns  YVAL  given a value  XVAL by interpolating   *
!   using a pair of arrays XARR(1:NARR), YARR(1:NARR), that tabulate   *
!   a  function.  ACCY  is the  desired  accuracy of the estimate: a   *
!   warning message is issued if this is not achieved.  A warning is   *
!   also issued when the routine is extrapolating.  Aitken's algori-   *
!   thm is used. See, for  instance, F B Hildebrand, Introduction to   *
!   Numerical Analysis,  2nd ed., McGraw-Hill, New York, NY, 1974.     *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: narr
      real(real64), intent(in) :: xval
      real(real64), intent(out) :: yval
      real(real64), intent(in) :: accy
      real(real64), intent(in) :: xarr(narr)
      real(real64), intent(in) :: yarr(narr)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: mxord = 11
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nrstlo, nrsthi, k, llo, lhi, llr, irow, locnxt, ilirok, ildiag&
         , ilothr, ibest, norder
      real(real64) :: dx(mxord), x(mxord), est(mxord)
      real(real64) :: poly((mxord*(mxord + 1))/2)
      real(real64) :: diff, difft, debeb, debe
      logical :: used(2*mxord + 2)
!-----------------------------------------------
!
!
!   MXORD is the maximum order of the interpolation
!
!
!   This function dispenses with the need for a two-dimensional
!   array for the Aitken lower triangle matrix
!
!
!   Determine the nearest two XARR entries bounding XVAL
!
      if (narr < 1) error stop "INTERP: empty table"
      norder = min(mxord,narr)
      nrstlo = 1
      if (xval < xarr(1)) then
         nrstlo = 1
         nrsthi = 1
         write (*, 300)
      else if (xval > xarr(narr)) then
         nrstlo = narr
         nrsthi = narr
         write (*, 300)
      else
         k = 0
    1    continue
         k = k + 1
         if (xarr(k) < xval) then
            nrstlo = k
            goto 1
         else
            nrsthi = k
         end if
      end if
      if (xval == xarr(nrsthi)) then
         yval = yarr(nrsthi)
         return
      end if
!
!   Clear relevant piece of use-indicator array
!
      llo = max(nrstlo - mxord,1)
      lhi = min(nrsthi + mxord,narr)
      llr = llo - 1
      used(llo-llr:lhi-llr) = .false.
!
!   Determine next-nearest XARR entry
!
      do irow = 1, norder
         llo = max(nrstlo - irow + 1,1)
         lhi = min(nrsthi + irow - 1,narr)
! Start from an actual unused sample, never from an undefined index.
         locnxt = llo
         do while (locnxt <= lhi)
            if (.not.used(locnxt-llr)) exit
            locnxt = locnxt + 1
         end do
         if (locnxt > lhi) error stop "INTERP: exhausted stencil"
         diff = xarr(locnxt) - xval
         do k = locnxt + 1, lhi
            if (used(k-llr)) cycle
            difft = xarr(k) - xval
            if (abs(difft) < abs(diff)) then
               diff = difft
               locnxt = k
            end if
         end do
         used(locnxt-llr) = .true.
         x(irow) = xarr(locnxt)
         dx(irow) = diff
!
!   Fill table for this row
!
         poly(iloc(irow,1)) = yarr(locnxt)
         do k = 2, irow
            ilirok = iloc(irow,k)
            ildiag = iloc(k - 1,k - 1)
            ilothr = iloc(irow,k - 1)
            poly(ilirok) = (poly(ildiag)*dx(irow)-poly(ilothr)*dx(k-1))/(x(&
               irow)-x(k-1))
         end do
!
!   Pick off the diagonal element
!
         ildiag = iloc(irow,irow)
         est(irow) = poly(ildiag)
!
      end do
!
!   Now the estimate vector is filled in, so obtain the
!   best estimate
!
      if (norder == 1) then
         yval = est(1)
         write (*,"(A)") "INTERP: One-point extrapolation; accuracy cannot be estimated."
         return
      end if
      debeb = relative_change(est(2),est(1))

      ibest = 2
      do irow = 3, norder
         debe = relative_change(est(irow),est(irow-1))

!ps 13/12/2017
!ps If NARR is small, the EST array sometimes contains the NaN values.
!ps In a consequence, DEBE becomes NaN, and in that case the condition
!ps in the CYCLE statement (NaN >= number) is then evaluated to "false".
!ps So, the CYCLE is not performed, what is incorrect.
!ps To avoid this, the CYCLE was replaced with the IF / THEN / END IF
!ps (the condition was reversed).
        !ps IF (DEBE >= DEBEB) CYCLE
        !ps DEBEB = DEBE
        !ps IBEST = IROW
        if (debe < debeb) then
            debeb = debe
            ibest = irow
        end if

      end do
      yval = est(ibest)
!
      if (debeb > accy) write (*, 301) debeb, accy
!
      return
!
  300 format("INTERP: Extrapolating, not interpolating.")
  301 format("INTERP: Accuracy of interpolation (",1p,1d10.3,") is",&
         " below input criterion (",1d10.3,").")
      return
      contains

! Keep the original relative criterion for nonzero estimates. Identical
! zero estimates have zero change; a new zero after a nonzero estimate has
! no finite relative accuracy estimate, so it must not silently pass.
      real(real64) function relative_change(current, previous) result(change)
      real(real64), intent(in) :: current, previous
      if (current == 0.0_real64) then
         change = huge(1.0_real64)
         if (previous == 0.0_real64) change = 0.0_real64
      else
         change = abs((current-previous)/current)
      end if
      end function relative_change

      integer function iloc (ind1, ind2)
      integer, intent(in) :: ind1
      integer, intent(in) :: ind2
      iloc = (ind1*(ind1 - 1))/2 + ind2
      return
      end function iloc
!
      end subroutine interp
end module interp_i
