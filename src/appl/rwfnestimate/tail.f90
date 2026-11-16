

!***********************************************************************
!                                                                      *
      subroutine tail(iorb, p, q, jp, mtp)
!                                                                      *
!   This subroutine begins the inward integration of the homogeneous   *
!   Dirac radial equation. With only minor modifications, the series   *
!   given by  J E Sienkiewicz  and W E Baylis, J Phys B: At Mol Phys   *
!   20 (1987) 5145-5156, p 5155, is used.                              *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use parameter_def, only: nnnp
      use def_c, only: accy, c, expmax
      use grid_c
      use orb_c
      use pote_c, only: yp
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: iorb
      integer , intent(in) :: jp
      integer , intent(out) :: mtp
      real(double) , intent(inout) :: p(nnnp)
      real(double) , intent(out) :: q(nnnp)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, nm4, m, loc
      real(double) :: eps, bige, bigebc, bgebc2, beta, fk, rjp, qem, offset, t&
         , zloc, gam2, zlbb, dnu, fkmzbb, sump, sumq, em, cm, emfact, ovltrm, &
         pterm, qterm, exptrm
!-----------------------------------------------
!
!
!   Initialize
!
      eps = accy*0.1d00
      bige = -e(iorb)
      bigebc = bige/c
      bgebc2 = bigebc/c
      beta = sqrt((-bige*(2.0d00 + bgebc2)))
      fk = dble(nak(iorb))
!
!   Find MTP
!
      i = jp
      nm4 = n - 4
      rjp = r(jp)
      qem = 0.25d00*expmax
    1 continue
      i = i + 1
      if (i <= nm4) then
         if (beta*(r(i)-rjp) > qem) then
            mtp = i
         else
            go to 1
         endif
      else
         write (*, 300)
         mtp = nm4
      endif
!
!   Compute offset for exponential function
!
      offset = beta*r(mtp)
!
!   Tabulate tail points
!
      do i = mtp, n
!
         t = 2.0d00*beta*r(i)
!
         zloc = yp(i)
         gam2 = fk**2 - (zloc/c)**2
         zlbb = zloc/beta
         dnu = zlbb*(1.0d00 + bgebc2)
         fkmzbb = fk - zlbb
!
         m = -1
         sump = 0.0d00
         sumq = 0.0d00
         m = m + 1
         em = dble(m)
         if (m == 0) then
            cm = 1.0d00
            emfact = 1.0d00
         else
            cm = cm*(gam2 - (dnu - em)**2)
            emfact = emfact*em
         endif
         ovltrm = cm*(t**(dnu - em)/emfact)
         pterm = ovltrm*(fkmzbb + em)*beta
         qterm = ovltrm*(fkmzbb - em)*bigebc
         sump = sump + pterm
         sumq = sumq + qterm
         do while(abs(pterm/sump)>=eps .or. abs(qterm/sumq)>=eps)
            m = m + 1
            em = dble(m)
            if (m == 0) then
               cm = 1.0d00
               emfact = 1.0d00
            else
               cm = cm*(gam2 - (dnu - em)**2)
               emfact = emfact*em
            endif
            ovltrm = cm*(t**(dnu - em)/emfact)
            pterm = ovltrm*(fkmzbb + em)*beta
            qterm = ovltrm*(fkmzbb - em)*bigebc
            sump = sump + pterm
            sumq = sumq + qterm
         end do
         exptrm = exp((-0.5d00*t) + offset)
         p(i) = sump*exptrm
         q(i) = sumq*exptrm
         if (p(i) /= 0.0d00) cycle
         loc = i + 1
         go to 4
      end do
      loc = n + 1
!
    4 continue
      p(loc:n) = 0.0d00
      q(loc:n) = 0.0d00
!
      return
!
  300 format('TAIL: Grid may be of insufficient extent')
      return
!
      end subroutine tail
