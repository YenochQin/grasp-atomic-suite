!***********************************************************************
!                                                                      *
      subroutine count(fr, mtpfr, nncff, sgn)
!                                                                      *
!   This subroutine counts the nodes in the radial function FR using   *
!   the criteria  given by C Froese Fischer, Comp Phys Rep, 3 (1986)   *
!   314-315 . The  function FR is assumed defined on the first MTPFR   *
!   points of the radial grid. The sign of the function at the first   *
!   oscillation is also determined.                                    *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:46:58   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------                       *
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: nnnp
      use coun_c
      use def_c, only: accy
      use grid_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: mtpfr
      integer, intent(out) :: nncff
      real(double), intent(out) :: sgn
      real(double), dimension(nnnp), intent(in) :: fr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, dimension(nnnp) :: lcext
      integer :: next, i, loc, nstps
      integer :: trace_status
      character(len=8) :: trace_env
      logical, save :: trace_count = .false., trace_init = .false.
      real(double) :: ext, emx, abfri, test, threse, ablcl
      if (.not. trace_init) then
         call get_environment_variable('GRASP_TRACE_COUNT', trace_env, status=trace_status)
         trace_count = trace_status == 0 .and. (trace_env(1:1) == '1' .or. trace_env(1:1) == 'y')
         trace_init = .true.
      end if
!-----------------------------------------------
!
!
!   (1) Find all extrema in FR
!   (2) Find the maximum amplitudes of FR
!
      next = 1
      ext = 0.0d00
      lcext(1) = 1
      emx = 0.0d00
      do i = 2, mtpfr
         abfri = abs(fr(i))
         test = abs(sign(1.0d00,fr(i))+sign(1.0d00,fr(i-1)))
         if (test <= accy) then
            next = next + 1
            lcext(next) = 0
            ext = 0.0d00
         endif
         if (abfri > ext) then
            ext = abfri
            lcext(next) = i
         endif
         if (abfri <= emx) cycle
         emx = abfri
      end do
!
!   Eliminate oscillations with amplitude less than THRESH times
!   the maximum
!
      loc = 0
      threse = thresh*emx
    4 continue
      loc = loc + 1
      if (loc <= next) then
         if (lcext(loc) == 0) then
            ablcl = 0.0d00
         else
            ablcl = abs(fr(lcext(loc)))
         endif
         if (ablcl < threse) then
            next = next - 1
            nstps = next - loc
            lcext(loc:nstps+loc) = lcext(loc+1:nstps+1+loc)
            loc = loc - 1
         endif
         go to 4
      endif
!
!   Count changes of sign using the remaining oscillations
!
      nncff = 0
      do i = 2, next
         test = abs(sign(1.0d00,fr(lcext(i)))+sign(1.0d00,fr(lcext(i-1))))
         if (test > accy) cycle
         nncff = nncff + 1
      end do
!
!   Determine the position of the first oscillation, and the
!   sign of the function at this location
!
      sgn = sign(1.0d00,fr(lcext(1)))
      if (trace_count .and. mtpfr > 300 .and. count_trace_rank == 0) then
         write (734,'(A,I0,A,I0,A,ES16.8,A,ES16.8,A,I0)') 'COUNT_TRACE orbital=', &
              count_context, ' mtp=', mtpfr, ' emx=', emx, ' thre=', threse, ' nodes=', nncff
         write (734,'(A,I0)') 'COUNT_EXTREMA orbital=', count_context
         do i = 1, next
            if (lcext(i) > 0) write (734,'(I0,A,I0,A,ES16.8,A,ES16.8)') &
                 i, ',', lcext(i), ',', r(lcext(i)), ',', abs(fr(lcext(i)))
         end do
      end if
!
      return
      end subroutine count
