!***********************************************************************
!                                                                      *
      subroutine intrpqi(pa, qa, ma, ra, j, dnorm)
!                                                                      *
!   This  subprogram  interpolates  the  arrays  PA(1:MA), QA(1:MA),   *
!   tabulated on grid RA(1:MA) into the COMMON arrays PF(1:MF(J),J),   *
!   QF(1:MF(J),J). (Aitken's  algorithm is used. See F B Hildebrand,   *
!   Introduction  to  Numerical  Analysis, 2nd ed., McGraw-Hill, New   *
!   York, NY, 1974.) The orbital is renormalized.                      *
!                                                                      *
!   SUBROUTINEs called: RINT.                                          *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_radial_points
      use debug_c
      use def_c,           only: accy
      use grid_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rintii_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ma
      integer, intent(in) :: j
      real(real64) , intent(out) :: dnorm
      real(real64) , intent(in) :: pa(*), qa(*), ra(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: mxord = 13
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, mfj, nrstlo, kount, irow, k, nrsthi, llo, lhi, locnxt, &
         ilirok, ildiag, ilothr, mfjp1
      real(real64) :: x(mxord), dx(mxord)
      real(real64) :: polyp((mxord*(mxord + 1))/2), polyq((mxord*(mxord + 1))/2)
      real(real64) :: rama, rn, xbar, pestl, qestl, diff, difft, dxkmn1, dxirow&
         , factor, pestt, qestt, dpbp, dqbq, dnfac
      logical :: set
      logical :: used(max_radial_points)
!-----------------------------------------------
!
!   MXORD is the maximum order of the interpolation
!
!
!   This function dispenses with the need for a two-dimensional
!   array
!
!   Initialization
!
      rama = ra(ma)
      rn = radial_grid(radial_point_count)
!
!   This is always true in GRASP
!
      pfii(1,j) = 0.0_real64
      qfii(1,j) = 0.0_real64
!
!   Checks
!
      if (rama > rn) then
         write (*, 300) rn, rama
         stop
      endif
!
!   Determine end of grid
!
      i = radial_point_count
    1 continue
      i = i - 1
      if (radial_grid(i) <= rama) then
         mfj = i
      else
         go to 1
      endif
      mfii(j) = mfj
!
!   Overall initialization for interpolation
!
      nrstlo = 0
      kount = 0
!
!   Perform interpolation
!
      do i = 2, mfj
!
!   Initialization for interpolation
!
         xbar = radial_grid(i)
         irow = 0
         pestl = 0.0_real64
         qestl = 0.0_real64
!
!   Determine the nearest two grid points bounding the present
!   grid point
!
    2    continue
         k = nrstlo + 1
         if (ra(k) < xbar) then
            nrstlo = k
            go to 2
         else
            nrsthi = k
         endif
!
!   Clear relevant piece of use-indicator array
!
         llo = max(nrstlo - mxord,1)
         lhi = min(nrsthi + mxord,ma)
         used(llo:lhi) = .false.
!
!   Determine next-nearest grid point
!
    4    continue
         irow = irow + 1
         llo = max(nrstlo - irow + 1,1)
         lhi = min(nrsthi + irow - 1,ma)
         set = .false.
         do k = llo, lhi
            if (used(k)) cycle
            if (.not.set) then
               diff = ra(k) - xbar
               locnxt = k
               set = .true.
            else
               difft = ra(k) - xbar
               if (abs(difft) < abs(diff)) then
                  diff = difft
                  locnxt = k
               endif
            endif
         end do
         used(locnxt) = .true.
         x(irow) = ra(locnxt)
         dx(irow) = diff
!
!   Fill table for this row
!
         do k = 1, irow
            ilirok = iloc(irow,k)
            if (k == 1) then
               polyp(ilirok) = pa(locnxt)
               polyq(ilirok) = qa(locnxt)
            else
               ildiag = iloc(k - 1,k - 1)
               ilothr = iloc(irow,k - 1)
               dxkmn1 = dx(k-1)
               dxirow = dx(irow)
               factor = 1.0_real64/(x(irow)-x(k-1))
               polyp(ilirok) = (polyp(ildiag)*dxirow-polyp(ilothr)*dxkmn1)*&
                  factor
               polyq(ilirok) = (polyq(ildiag)*dxirow-polyq(ilothr)*dxkmn1)*&
                  factor
            endif
         end do
!
!   Check for convergence
!
         ildiag = iloc(irow,irow)
         pestt = polyp(ildiag)
         qestt = polyq(ildiag)
         if (pestt==0.0_real64 .or. qestt==0.0_real64) then
            if (irow < mxord) then
               go to 4
            else
               pfii(i,j) = pestt
               qfii(i,j) = qestt
            endif
         else
            dpbp = abs((pestt - pestl)/pestt)
            dqbq = abs((qestt - qestl)/qestt)
            if (dqbq<accy .and. dpbp<accy) then
               pfii(i,j) = pestt
               qfii(i,j) = qestt
            else
               pestl = pestt
               qestl = qestt
               if (irow < mxord) then
                  go to 4
               else
                  pfii(i,j) = pestt
                  qfii(i,j) = qestt
                  kount = kount + 1
               endif
            endif
         endif
!
      end do
!
!   Ensure that all points of the array are defined by setting the
!   tail to zero
!
      mfjp1 = mfj + 1
      pfii(mfjp1:radial_point_count,j) = 0.0_real64
      qfii(mfjp1:radial_point_count,j) = 0.0_real64
!
      if (ldbpr(3) .and. kount>0) write (99, 301) accy, kount, mfj
!
!   Normalization
!
      dnorm = rintii(j,j,0)
!ww      WRITE(92,*) 'DNORM',DNORM
      dnfac = 1.0_real64/dsqrt(dnorm)
      pfii(:mfj,j) = pfii(:mfj,j)*dnfac
      qfii(:mfj,j) = qfii(:mfj,j)*dnfac
!
      return
!
  300 format(/,'INTRPQ: Grid of insufficient extent:'/,&
         ' Present grid has R(N) = ',1p,1d19.12,' Bohr radii'/,&
         '          Require R(N) = ',1d19.12,' Bohr radii')
  301 format(/,'INTRPQ: Interpolation procedure not converged to',1p,1d19.12,&
         ' for ',1i3,' of ',1i3,' tabulation points')
      return
      contains


      integer function iloc (ind1, ind2)
      integer, intent(in) :: ind1
      integer, intent(in) :: ind2
      iloc = (ind1*(ind1 - 1))/2 + ind2
      return
      end function iloc
!
      end subroutine intrpqi
