!***********************************************************************
!                                                                      *
      subroutine start(iorb, itype, p0, p, q0, q)
!                                                                      *
!   This subroutine sets up  P(1:6), Q(1:6),  required  to start the   *
!   integration for programs  OUT  and  SBSTEP .                       *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!      IORB : (Input) Index of the orbital                             *
!      ITYPE: (Input) 1 = homogeneous equation; 2 = inhomogeneous      *
!             equation; 3 = variational equation                       *
!      P0   : (Input) Slope parameter                                  *
!      P    : (Output) P(1:6) are tabulated by this program            *
!      Q0   : (Output) First term in the series expansion of Q         *
!      Q    : (Output) Q(1:6) are tabulated by this program            *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:48   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: nnnp
      use cnc_c
      use def_c,           only: c, z, accy
      use grid_c
      use npar_c
      use orb_c
      use pote_c,          only: yp, xp, xq
      use scf_c,           only: ndcof, nda, da
      use wave_c,          only: pz,pf,qf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: iorb
      integer, intent(in) :: itype
      real(double), intent(in) :: p0
      real(double), intent(out) :: q0
      real(double), dimension(nnnp), intent(inout) :: p
      real(double), dimension(nnnp), intent(inout) :: q
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mxiter=36, kiorb, i, jorb, niter, j
      real(double), dimension(6) :: rdp, rdq, rsep, rseq
      real(double), dimension(2:6) :: spest, sqest
      real(double) :: obc, zonc, giorb, fkiorb, omgi, omgmk, omgpk, rsep1, &
         rseq1, pzero, p1, q1, sump, sumq, factor, csq, twocsq, energy, enefac&
         , ri, rpi, rirpi, ypirpi, difmaw, difmax, coefij, pi, qi, rj
!-----------------------------------------------
!
!
!   Initialization
!
      obc = 1.0d00/c
      zonc = z*obc
      giorb = gama(iorb)
      kiorb = nak(iorb)
      fkiorb = dble(kiorb)
      omgi = 1.0d00 - giorb
!
      omgmk = omgi - fkiorb
      omgpk = omgi + fkiorb
!
!   Determine P(1), Q(1): THESE STORE  R**(-GAMMA)*(P(1),Q(1));
!   set up  RSEP  and  RSEQ , the inhomogeneous terms
!
      if (itype==1 .or. itype==2) then
         if (nparm == 0) then
            p(1) = p0
            if (kiorb < 0) then
               q(1) = -p0*zonc/(giorb - fkiorb)
            else
               q(1) = p0*(giorb + fkiorb)/zonc
            endif
         else
            if (kiorb < 0) then
               p(1) = p0
               q(1) = 0.0d00
            else
               p(1) = 0.0d00
               q(1) = p0*(giorb + fkiorb)/zonc
            endif
         endif
         if (itype == 1) then
            rsep = 0.0d00
            rseq = 0.0d00
         else
            rsep1 = 0.0d00
            rseq1 = 0.0d00
            do i = 1, ndcof
               jorb = nda(i)
               pzero = pz(jorb)
               if (nparm == 0) then
                  p1 = pzero
                  if (kiorb < 0) then
                     q1 = -pzero*zonc/(giorb - fkiorb)
                  else
                     q1 = pzero*(giorb + fkiorb)/zonc
                  endif
                  sump = zonc*q1
                  sumq = -zonc*p1
               else
                  if (kiorb < 0) then
                     p1 = pzero
                     q1 = 0.0d00
                  else
                     p1 = 0.0d00
                     q1 = pzero*(giorb + fkiorb)/zonc
                  endif
                  sump = 0.0d00
                  sumq = 0.0d00
               endif
               factor = da(i)
               rsep1 = rsep1 + factor*(sump + omgmk*p1)
               rseq1 = rseq1 + factor*(sumq + omgpk*q1)
            end do
            factor = rp(1)
            rsep(1) = factor*rsep1
            rseq(1) = factor*rseq1
            do i = 2, 6
               factor = -rp(i)*r(i)**(-giorb)
               rsep(i) = factor*xp(i)
               rseq(i) = factor*xq(i)
            end do
         endif
      else if (itype == 3) then
         p(1) = 0.0d00
         q(1) = 0.0d00
         rsep(1) = 0.0d00
         rseq(1) = 0.0d00
         do i = 2, 6
            factor = obc*rp(i)*r(i)**omgi
            rsep(i) = -factor*qf(i,iorb)
            rseq(i) = factor*pf(i,iorb)
         end do
      endif
      q0 = q(1)
!
!   Set up  RDP  and  RDQ
!
      csq = c*c
      twocsq = csq + csq
      energy = e(iorb)
      enefac = twocsq - energy
      do i = 1, 6
         ri = r(i)
         rpi = rp(i)
         rirpi = ri*rpi
         ypirpi = yp(i)*rpi
         rdp(i) = -obc*(enefac*rirpi + ypirpi)
         rdq(i) = -obc*(energy*rirpi - ypirpi)
      end do
!
!   Determine  P(2:6) , Q(2:6)
!
!   Initilizations for the iterations
!
      niter = 0
      p1 = p(1)
      q1 = q(1)
      difmaw = max(abs(p1),abs(q1))
!
      p(2:6) = p1
      q(2:6) = q1
!
!   This is the largest factor by which any result will be
!   multiplied
!
      factor = r(6)**giorb
!
!   Now iterate
!
    7 continue
      niter = niter + 1
      difmax = 0.0d00
      do j = 2, 6
         sump = sum(cnc6c(:,j)*(omgmk*rp(:6)*p(:6)-rdp*q(:6)+rsep))
         sumq = sum(cnc6c(:,j)*(omgpk*rp(:6)*q(:6)-rdq*p(:6)+rseq))
         rj = r(j)
         sump = sump/rj
         sumq = sumq/rj
         spest(j) = sump
         sqest(j) = sumq
         difmax = max(difmax,abs(sump - p(j)))
         difmax = max(difmax,abs(sumq - q(j)))
      end do
!zou  IF (DIFMAX .LT. DIFMAW) THEN
      p(2:6) = spest(:6)
      q(2:6) = sqest(:6)
      difmaw = difmax
      difmax = difmax*factor
      if (difmax > accy) then
         if (niter < mxiter) then
            go to 7
         endif
      endif
!     ELSE
!        DIFMAX = DIFMAX*FACTOR
!        IF (DIFMAX .GT. ACCY) THEN
!           WRITE (*,300) NP(IORB),NH(IORB),DIFMAX,NITER,ACCY
!        ENDIF
!zou  ENDIF
! not convergent, using the initial P,Q
      if (difmax > accy) then
         p(2:6) = p1
         q(2:6) = q1
      endif
!zou
!
!   All done
!
!   This is always true in GRASP2
!
      p(1) = 0.0d00
      q(1) = 0.0d00
!
      do i = 2, 6
         factor = r(i)**giorb
         p(i) = factor*p(i)
         q(i) = factor*q(i)
      end do
!
      return
!
  300 format('START: ',1i2,1a2,' subshell: accuracy ',1p,1d8.1,/,&
         ' attained after ',1i2,' iterations; this fails the'/,&
         ' accuracy criterion ',d8.1,'.')
      return
!
      end subroutine start
