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
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_orbitals
      use decide_c, only: double, lfordr, lnms, lse, lsms, ltrans, lvp
      use eigv_c
      use foparm_c,        only: iccut
      use iounit_c,        only: stderr_unit
      use orb_c,           only: csf_count, orbital_count
      use prnt_c,          only: nvec
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_i, only: convrt
      use itjpo_i, only: itjpo
      use matelt_i, only: matelt
      use oneparticlejj_i, only: oneparticlejj
      use rint_i, only: rint
      use rinthf_i, only: rinthf
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(out) :: gjc_diag(nvec), dgjc_diag(nvec)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(real64), parameter :: cutoff = 1.0e-10_real64
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, ipt, ic, lcnum, ir, itjpoc, itjpor, idiff, ia, ib, &
         k, loc
      real(real64) :: tshell(max_orbitals)
      real(real64) :: rintgj(max_orbitals, max_orbitals), rintdgj(max_orbitals, max_orbitals), &
         gjmelt(max_orbitals, max_orbitals), dgjmelt(max_orbitals, max_orbitals)
      real(real64) :: apart, gjpart, dgjpart, elemntgj, elemntdgj, contrgj, &
         contrdgj
      character :: cnum*11
!-----------------------------------------------
!
      gjc_diag(:nvec) = 0.0_real64
      dgjc_diag(:nvec) = 0.0_real64
!
!   Precompute the radial integrals and angular factors required by
!   the two one-body operators contributing to the Landé factor.
!
      do i = 1, orbital_count
         do j = 1, orbital_count
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
      do ic = 1, csf_count
         if (mod(ic,100) == 0) then
            call convrt(ic, cnum, lcnum)
            write (stderr_unit, *) "Column "//cnum(1:lcnum)//" complete;"
         end if
!
         do ir = 1, csf_count
            if (lfordr .and. ic>iccut .and. ic/=ir) cycle
!
            itjpoc = itjpo(ic)
            itjpor = itjpo(ir)
            idiff = itjpoc - itjpor
            if (.not.(idiff == 0 .and. ir >= ic)) cycle
!
            elemntgj = 0.0_real64
            elemntdgj = 0.0_real64
!
            call oneparticlejj(1, ipt, ic, ir, ia, ib, tshell)
            if (ia == 0) cycle
!
            if (ia == ib) then
               do i = 1, orbital_count
                  if (abs(tshell(i)) <= cutoff) cycle
                  elemntgj = elemntgj + gjmelt(i,i)*rintgj(i,i)*tshell(i)
                  elemntdgj = elemntdgj + dgjmelt(i,i)*rintdgj(i,i)*tshell(i)
               end do
            else
               if (abs(tshell(1)) <= cutoff) cycle
               elemntgj = elemntgj + gjmelt(ia,ib)*rintgj(ia,ib)*tshell(1)
               elemntdgj = elemntdgj + dgjmelt(ia,ib)*rintdgj(ia,ib)*tshell(1)
            end if
!
            if (abs(elemntgj) <= cutoff .and. abs(elemntdgj) <= cutoff) cycle
!
            do k = 1, nvec
               loc = (k - 1)*csf_count
               if (ir /= ic) then
                  contrgj = elemntgj*(csf_mixing_coefficients(ic + loc)*csf_mixing_coefficients(ir + loc) + &
                     & csf_mixing_coefficients(ir &
                     + loc)*csf_mixing_coefficients(ic + loc))
                  contrdgj = elemntdgj*(csf_mixing_coefficients(ic + loc)*csf_mixing_coefficients(ir + loc) + &
                     & csf_mixing_coefficients(ir&
                      + loc)*csf_mixing_coefficients(ic + loc))
               else
                  contrgj = elemntgj*csf_mixing_coefficients(ic + loc)*csf_mixing_coefficients(ir + loc)
                  contrdgj = elemntdgj*csf_mixing_coefficients(ic + loc)*csf_mixing_coefficients(ir + loc)
               end if
               gjc_diag(k) = gjc_diag(k) + contrgj
               dgjc_diag(k) = dgjc_diag(k) + contrdgj
            end do
         end do
      end do
!
      return
      end subroutine compute_gj_only
