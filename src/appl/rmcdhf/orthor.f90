!***********************************************************************
!                                                                      *
      subroutine orthor(j, inv)
!
!                                                                      *
!   This routine Schmidt orthogonalizes orbital  J  to  all orbitals   *
!   which  have  better  self-consistency.  Note that fixed orbitals   *
!   have the best self-consistency.                                    *
!                                                                      *
!   Call(s) to: [LIB92]: COUNT, RINT.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
!                                                                      *
!***********************************************************************
! INV not initialized here but used as
! INV = INV + 1
! A bug regarding IORDER fixed and new schemes introduced. See comments
! Anyway this routine was not used anywhere.
!XHH 1997.02.21
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:54:29   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: nnnp
      use def_c
      use grid_c
      use ovl_c
      use orb_c
      use orba_c
      use invt_c
      use scf_c
      use iounit_c
      use fixd_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rint_i
!      USE count_I
      use count_i, only: count
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
      integer , intent(inout) :: inv
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kraw, k, mtp, i, nncff, jgg
      real(double) :: eps, ovrlap, dnorm, factor, sgn
      logical :: check, changed
!-----------------------------------------------
!
!XHH Added common /orba/ and /iounit/

!
!ww      EPS = ACCY*0.1D 00
      eps = accy*0.01d00
!
      check = .not.noinvt(j)
!
!XHH
! Bug fixed. Here J is actually IORDER(J_raw), thus K should be
! treated the same way.
      changed = .false.
!      DO 4 K = 1,NW
      do kraw = 1, nw
         k = iorder(kraw)
!         write(istde,*) '***',kraw,k,np(k),nh(k),scnsty(k),'***'
!
!XHH orbitals with higher self-consistency are considered
!         IF ( (K .NE. J) .AND.
!     :        (NAK(K) .EQ. NAK(J)) .AND.
!     :        (SCNSTY(K) .LT. SCNSTY(J)) ) THEN
!XHH All orbitals are considered
!         IF ( (K .NE. J) .AND.
!     :        (NAK(K) .EQ. NAK(J)) ) THEN
!XHH orbitals before the current and unchanged ones are considered
         if (.not.(nak(k)==nak(j) .and. (k<j .or. lfix(k)))) cycle

         changed = .true.
!
!   Compute overlap
!
         ovrlap = rint(j,k,0)
!
!   Schmidt orthogonalise
!
         pz(j) = pz(j) - ovrlap*pz(k)
         mtp = max(mf(j),mf(k))
         pf(:mtp,j) = pf(:mtp,j) - ovrlap*pf(:mtp,k)
         qf(:mtp,j) = qf(:mtp,j) - ovrlap*qf(:mtp,k)
      end do
!XHH
      if (changed) then
!
!   Normalise
!
         mf(j) = mtp
         dnorm = rint(j,j,0)
         factor = 1.0d00/sqrt(dnorm)
!
!   Determine if inversion is necessary
!
         if (check) then
            call count (pf(:nnnp,j), mtp, nncff, sgn)
            if (sgn < 0.0d00) then
               inv = inv + 1
               factor = -factor
            endif
         endif
!
!   Perform normalization and/or inversion
!
         pz(j) = factor*pz(j)
         pf(2:mtp,j) = factor*pf(2:mtp,j)
         qf(2:mtp,j) = factor*qf(2:mtp,j)
!
!   Find new MF(J)
!
         mtp = mtp + 1
    3    continue
         mtp = mtp - 1
         if (abs(pf(mtp,j)) < eps) then
            pf(mtp,j) = 0.0d00
            qf(mtp,j) = 0.0d00
            go to 3
         else
            mf(j) = mtp
         endif
!
!XHH two lines moved forward and added an ENDIF
!         ENDIF
!    4 CONTINUE
      endif
!
      return
      end subroutine orthor
