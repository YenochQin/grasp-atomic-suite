!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function vinti (j, k)
!                                                                      *
!   The value of this  function is the one-electron integral V (J,K)   *
!   for  orbitals  J, K. The analytical expression for this quantity   *
!   is given as  eq (3.23) in  F A  Parpia, M. Tong and C F Fischer,   *
!   to appear.                                                         *
!                                                                      *
!   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
!                                                                      *
!   Written by M Tong and F A Parpia,     Last revision: 15 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:07:11   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use debug_c
      use grid_c
      use orb_c
      use tatb_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dpbdt_i
      use quad_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
      integer  :: k
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, kpj, kpk, ifact1, ifact2
      real(double) :: piece1, fact1, fact2, piece2
!-----------------------------------------------
!
      mtp = max(mf(j),mf(k))
!
!   Piece involving derivatives
!
      call dpbdt (k)
      ta(1) = 0.0d00
      do i = 2, mtp
         ta(i) = pf(i,j)*ta(i) + qf(i,j)*tb(i)
      end do
      call quad (piece1)
      piece1 = piece1/h
!
!   Pieces not involving derivatives
!
      kpj = nak(j)
      kpk = nak(k)
      ifact1 = kpj*(kpj + 1) - kpk*(kpk + 1)
      fact1 = 0.5d00*dble(ifact1)
      ifact2 = (-kpj*((-kpj) + 1)) + kpk*((-kpk) + 1)
      fact2 = 0.5d00*dble(ifact2)
      ta(1) = 0.0d00
      do i = 2, mtp
         ta(i) = rpor(i)*(fact1*pf(i,j)*pf(i,k) + fact2*qf(i,j)*qf(i,k))
      end do
      call quad (piece2)
!
      vinti = piece1 - piece2
!
!   Debug printout
!
      if (ldbpr(6)) write (99, 300) np(j), nh(j), np(k), nh(k), vinti
!
      return
!
  300 format(/,'VINTI: V (',1i2,1a2,',',1i2,1a2,') = ',1p,d19.12)
      return
!
      end function vinti
