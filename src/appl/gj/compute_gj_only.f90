!***********************************************************************
!                                                                      *
      subroutine compute_gj_only(gjc_diag, dgjc_diag)
!                                                                      *
!   This routine evaluates only the diagonal ASF contributions needed  *
!   for g_J, delta g_J, and total g_J.                                 *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: nnnw
      use decide_c
      use eigv_c
      use foparm_c,        only: iccut
      use iounit_c,        only: istde
      use orb_c,           only: ncf, nw
      use prnt_c,          only: nvec
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_i
      use itjpo_i
      use matelt_i
      use oneparticlejj_i
      use rint_i
      use rinthf_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double), intent(out) :: gjc_diag(nvec), dgjc_diag(nvec)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: cutoff = 1.0d-10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, ipt, ic, lcnum, ir, itjpoc, itjpor, idiff, ia, ib, &
         k, loc
      real(double), dimension(nnnw) :: tshell
      real(double), dimension(nnnw,nnnw) :: rintgj, rintdgj, gjmelt, dgjmelt
      real(double) :: apart, gjpart, dgjpart, elemntgj, elemntdgj, contrgj, &
         contrdgj
      character :: cnum*11
!-----------------------------------------------
!
      gjc_diag(:nvec) = 0.0d00
      dgjc_diag(:nvec) = 0.0d00
!
!   Precompute the radial integrals and angular factors required by
!   the two one-body operators contributing to the Landé factor.
!
      do i = 1, nw
         do j = 1, nw
            rintgj(i,j) = rinthf(i,j,1)
            rintdgj(i,j) = rint(i,j,0)
            call matelt(i, 1, j, apart, gjpart, dgjpart)
            gjmelt(i,j) = gjpart
            dgjmelt(i,j) = dgjpart
         end do
      end do
!
      ipt = 1
!
!   Sweep over the CSF pairs. Only the diagonal-in-J branch from HFSGG
!   is required, and only ASF diagonal matrix elements are accumulated.
!
      do ic = 1, ncf
         if (mod(ic,100) == 0) then
            call convrt(ic, cnum, lcnum)
            write (istde, *) 'Column '//cnum(1:lcnum)//' complete;'
         endif
!
         do ir = 1, ncf
            if (lfordr .and. ic>iccut .and. ic/=ir) cycle
!
            itjpoc = itjpo(ic)
            itjpor = itjpo(ir)
            idiff = itjpoc - itjpor
            if (.not.(idiff == 0 .and. ir >= ic)) cycle
!
            elemntgj = 0.0d00
            elemntdgj = 0.0d00
!
            call oneparticlejj(1, ipt, ic, ir, ia, ib, tshell)
            if (ia == 0) cycle
!
            if (ia == ib) then
               do i = 1, nw
                  if (abs(tshell(i)) <= cutoff) cycle
                  elemntgj = elemntgj + gjmelt(i,i)*rintgj(i,i)*tshell(i)
                  elemntdgj = elemntdgj + dgjmelt(i,i)*rintdgj(i,i)*tshell(i)
               end do
            else
               if (abs(tshell(1)) <= cutoff) cycle
               elemntgj = elemntgj + gjmelt(ia,ib)*rintgj(ia,ib)*tshell(1)
               elemntdgj = elemntdgj + dgjmelt(ia,ib)*rintdgj(ia,ib)*tshell(1)
            endif
!
            if (abs(elemntgj) <= cutoff .and. abs(elemntdgj) <= cutoff) cycle
!
            do k = 1, nvec
               loc = (k - 1)*ncf
               if (ir /= ic) then
                  contrgj = elemntgj*(evec(ic + loc)*evec(ir + loc) + evec(ir &
                     + loc)*evec(ic + loc))
                  contrdgj = elemntdgj*(evec(ic + loc)*evec(ir + loc) + evec(ir&
                      + loc)*evec(ic + loc))
               else
                  contrgj = elemntgj*evec(ic + loc)*evec(ir + loc)
                  contrdgj = elemntdgj*evec(ic + loc)*evec(ir + loc)
               endif
               gjc_diag(k) = gjc_diag(k) + contrgj
               dgjc_diag(k) = dgjc_diag(k) + contrdgj
            end do
         end do
      end do
!
      return
      end subroutine compute_gj_only
