!***********************************************************************
!                                                                      *
      subroutine dampor(j, inv, odampj)
!                                                                      *
!   This subroutine damps the orbital wave function with index J. it   *
!   also  stores  the  previous  determination  of this orbital.       *
!                                                                      *
!   Call(s) to: [LIB92]: COUNT, RINT.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:22:36   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: nnnp
      use damp_c
      use def_c
      use grid_c
      use int_c
      use invt_c
      use scf_c
      use wave_c
      use mpi_c
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
      integer, intent(inout) :: inv
      real(double), intent(in) :: odampj
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mtpo, mtpn, mtp, i, nncff, mfj
      real(double) :: eps, factor, poldi, qoldi, dnorm, dnfac, sgn
      logical :: check
!-----------------------------------------------
!
!
!   Initialization
!
!ww      EPS = 0.1D 00*ACCY
      eps = 0.01d00*accy
      check = .not.noinvt(j)
!
!   Damp orbital J using the damping factor  ABS (ODAMP(J)):  ODAMP(J)
!   is restricted to the open interval (-1,1) by  DATSCF ; the meaning
!   of a  negative  ODAMP(J)  is that  ABS (ODAMP(J))  is the constant
!   damping factor;  a positive  ODAMP(J)  can be reset by  subroutine
!   dampck
!
!   Store previous determination of this orbital in arrays P, Q
!
!XHH odampj goes to the argument
!      ODAMPJ = ABS (ODAMP(J))
!
      if (odampj > eps) then
!
         factor = 1.0d00 - odampj
!
         pz(j) = factor*p0 + odampj*pz(j)
!
         mtpo = mf(j)
         mtpn = mtp0
         mtp0 = mtpo
!
         mtp = max(mtpn,mtpo)
         do i = 1, mtp
            poldi = pf(i,j)
            pf(i,j) = factor*p(i) + odampj*pf(i,j)
            p(i) = poldi
            qoldi = qf(i,j)
            qf(i,j) = factor*q(i) + odampj*qf(i,j)
            q(i) = qoldi
         end do
!
!   Compute normalization factor
!
         mf(j) = mtp
         dnorm = rint(j,j,0)
         dnfac = 1.0d00/sqrt(dnorm)
!
!   Determine if inversion is necessary
!
         if (check) then
            call count (pf(:nnnp,j), mtp, nncff, sgn)
            if (sgn < 0.0d00) then
               inv = inv + 1
               dnfac = -dnfac
            endif
         endif
!
         pz(j) = pz(j)*dnfac
         pf(:mtp,j) = dnfac*pf(:mtp,j)
         qf(:mtp,j) = dnfac*qf(:mtp,j)
!
!   Find new MF(J)
!
         mfj = mtp + 1
    3    continue
         mfj = mfj - 1
         if (abs(pf(mfj,j)) < eps) then
            pf(mfj,j) = 0.0d00
            qf(mfj,j) = 0.0d00
            go to 3
         else
            mf(j) = mfj
         endif
!
      else
!
         pz(j) = p0
!
         mtpo = mf(j)
         mtpn = mtp0
         mtp0 = mtpo
!
         mtp = max(mtpn,mtpo)
         do i = 1, mtp
            poldi = pf(i,j)
            pf(i,j) = p(i)
            p(i) = poldi
            qoldi = qf(i,j)
            qf(i,j) = q(i)
            q(i) = qoldi
         end do
!
         mf(j) = mtp
!
      endif
!
      return
      end subroutine dampor
