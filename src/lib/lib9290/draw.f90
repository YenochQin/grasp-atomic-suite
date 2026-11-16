!***********************************************************************
!                                                                      *
      subroutine draw(p, sp, q, sq, mf)
!                                                                      *
!   This  subroutine  generates  a  printer plot. P and Q are radial   *
!   functions with the maximum tabulation point MF. SP is the factor   *
!   by which P is to be scaled, SQ is the factor by which Q is to be   *
!   scaled.                                                            *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 10 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_radial_points
      use grid_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: mf
      real(real64), intent(in) :: sp
      real(real64), intent(in) :: sq
      real(real64), intent(in) :: p(max_radial_points)
      real(real64), intent(in) :: q(max_radial_points)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ioffst, ixlocq, ixlocp, iyloc
      real(real64) :: sixty, fift4, ohtht, dmx, dmn, spi, sqi, rmx, xscal, &
         yscal, sqx, spx
      logical :: first
      character :: cblank*132, cdash*132
!cjb
!     CHARACTER , DIMENSION(60) :: CPLOT*132
      character  :: cplot(mf)*132
!cjb
      character :: cp, cq
!-----------------------------------------------
      data sixty/ 60.0_real64/
      data fift4/ 54.0_real64/
      data ohtht/ 131.0_real64/
!
      data first/ .true./


!
!   This initialization is carried out once per run
!
      if (first) then
         cblank(1:1) = '|'
         cdash(1:1) = '-'
         do i = 2, 131
            cblank(i:i) = ' '
            cdash(i:i) = '-'
         end do
         cblank(132:132) = '|'
         cdash(132:132) = '-'
         first = .false.
      endif
!
!   Initialization
!
      if (sq == 0.0_real64) then
         cp = 'X'
         cq = '-'
      else
         cp = 'P'
         cq = 'Q'
      endif
!
!   Determine the range of amplitude
!
      dmx = max(sp*p(1),sq*q(1))
      dmn = min(sp*p(1),sq*q(1))
      do i = 2, mf
         spi = sp*p(i)
         sqi = sq*q(i)
         dmx = max(max(dmx,spi),sqi)
         dmn = min(min(dmn,spi),sqi)
      end do
!
!   Determine the radial extent of the function
!
      rmx = radial_grid(mf)
!
!   Determine the scale factors
!
      if (dmx==0.0_real64 .and. dmn==0.0_real64) then
         write (99, 300)
         return
      else
         xscal = fift4/abs(dmx - dmn)
      endif
      yscal = ohtht/dble(radial_point_count)
!
!   Locate x = 0 if this is in the range
!
      if (dmx>0.0_real64 .and. dmn<0.0_real64) then
         ioffst = 4 + xscal*abs(dmn)
      else
         ioffst = 1
      endif
!
!   Initialize the array CPLOT
!
      do i = 1, 60
         if (i==1 .or. i==ioffst .or. i==60) then
            cplot(i) = cdash
         else
            cplot(i) = cblank
         endif
      end do
!
!   Generate the plot
!
!   Note that 'P' is the dominant character
!
      sqx = sq*xscal
      spx = sp*xscal
      do i = 1, mf
         ixlocq = nint(sqx*q(i)) + ioffst
         ixlocp = nint(spx*p(i)) + ioffst
         iyloc = max(1,nint(dble(i)*yscal))
         cplot(ixlocq)(iyloc:iyloc) = cq
         cplot(ixlocp)(iyloc:iyloc) = cp
      end do
!
!   Print plot
!
      write (99, 301)
      do i = 60, 1, -1
         write (99, 302) cplot(i)
      end do
!
!   Plot-information line
!
      write (99, 303) rmx, dmx
!
      return
!
  300 format(' No plot: function is identically zero')
  301 format('1')
  302 format(1x,1a132)
  303 format(/,50x,1p,' r (max) = ',1d10.3,'Bohr radii,',&
         ' Maximum of functions = ',1d10.3)
      return
!
      end subroutine draw
