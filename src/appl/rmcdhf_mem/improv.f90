!***********************************************************************
      subroutine improv(eol, j, lsort, dampmx)
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use corre_c
      use damp_c
      use def_c
      use grid_c
      use int_c
      use mpi_s
      use orb_c
      use orthct_c
      use scf_c
      use tatb_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use cofpot_i
      use defcor_i
      use solve_i
      use orthsc_i
      use matrix_i
      use newco_i
      use setlag_i
      use quad_i
      use consis_i
      use dampck_i
      use dampor_i
      use orthy_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
      real(double) , intent(inout) :: dampmx
      logical  :: eol, lsort
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: p2 = 2.0d-01
      real(double), parameter :: p005 = 5.0d-03
      real(double), parameter :: p0001 = 1.0d-04
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ipr, npts, inv, jp, nnp, i, nwww
      real(double) :: ed1, gamaj, ed2, wtaev, dnorm, dnfac, del1, del2, odampj
      logical :: fail, first
!-----------------------------------------------
!
!   C Froese Fischer's IPR and ED1 parameter
!
      data ipr/ 0/
      data ed1/ 0.d0/
      data first/ .false./
!
!
!-----------------------------------------------------------------------
      gamaj = gama(j)
!
!   C Froese Fischer's parameters IPR, ED1, ED2 are set and
!   used in this routine and in DAMPCK
!
    1 continue
      ed2 = e(j)
      ed1 = ped(j)
!
!   Set up the exchange potential and arrays XU, XV as appropriate
!
!   Set coefficients for YPOT, XPOT, DACON
!   Compute direct potential, exchange potential
!   Add in Lagrange-multiplier contribution
!   Add in derivative-terms contribution
!
      npts = n
      call cofpot (eol, j, npts)
!
!   Calculate deferred corrections
!
      call defcor (j)
!
!   Solve the Dirac equation
!
      inv = 0
      call solve (j, fail, inv, jp, nnp)
!
!   Upon failure issue message; take corrective action if possible
!
      if (fail) then
         if (myid == 0) write (*, 300) np(j), nh(j), method(j)
         if (method(j) /= 2) then
            method(j) = 2
!XHH orthsc does not have any argument
!    Orbital J [PF() and QF()]is not updated, why redo orthogonalization
            call orthsc
!CFF        ... avoid rediagonalization
!            IF (EOL) THEN
!               CALL MATRIX
!               CALL NEWCO (WTAEV)
!            ENDIF
            call setlag (eol)
            go to 1
         else
            if (myid == 0) write (*, 301)
            !CALL TIMER (0)
            stop
         endif
      endif
!
!   Compute norm of radial function
!
      ta(1) = 0.d0
      ta(2:mtp0) = (p(2:mtp0)**2+q(2:mtp0)**2)*rp(2:mtp0)
      mtp = mtp0

      call quad (dnorm)

!   Determine self-consistency [multiplied by SQRT(UCF(J))]

      call consis (j)
!
!   Normalize
!
      dnfac = 1.d0/sqrt(dnorm)
      p0 = p0*dnfac
      p(:mtp0) = p(:mtp0)*dnfac
      q(:mtp0) = q(:mtp0)*dnfac
!
!   Check if different method should be used or if improvement
!   count should be reduced
!
      del1 = abs(1.d0 - ed2/e(j))
      if (method(j) == 1) then
         del2 = max(abs(1.d0 - sqrt(dnorm)),abs(dnfac - 1.d0))
         if (del1<p005 .and. del2>p2) then
            method(j) = 2
            go to 1
         endif
      else
         if (del1<p0001 .and. nsic>1) nsic = nsic - 1
      endif
!
!   Damp the orbital --- if not converged
!
      if (scnsty(j) > accy) then
         call dampck (ipr, j, ed1, ed2)
         odampj = abs(odamp(j))
      else
         odampj = 0.d0                           ! take the whole new orbital
      endif
      call dampor (j, inv, odampj)

!   Orthogonalize all orbitals of the same kappa in the order
!   fixed, spectroscopic, correlation orbitals. The order of
!   orbitals in the latter two classes are sorted according
!   to their self-consistency and energy.

      if (orthst) then
         !CALL orthor (J, inv)
         nwww = nw
         call orthy (nwww, j, lsort)
      endif
!
!   Print details of iteration
!
      if (myid == 0)                                                &
         write (*, 302) np(j),nh(j),e(j),method(j),pz(j),scnsty(j), &
!cjb DNORM-1 -> SQRT(DNORM)-1
!cjb                    DNORM - 1, ODAMPJ, JP, MF(J), INV, NNP
                        sqrt(dnorm)-1,odampj,jp,mf(j),inv,nnp
      dampmx = max(dampmx,abs(odampj))

  300 format(/,' Failure; equation for orbital ',1i2,1a2,&
         ' could not be solved using method ',1i1)
  301 format(/,/,' ****** Error in SUBROUTINE IMPROV ******'/,&
         ' Convergence not obtained'/)
  302 format (1x,1i2,1a2,1p,1d16.7,1x,1i2,d11.3,1d10.2,1d10.2,&
              0p,f6.3,1x,1i5,1x,1i5,1x,1i2,1x,1i2)
      return
      end subroutine improv
