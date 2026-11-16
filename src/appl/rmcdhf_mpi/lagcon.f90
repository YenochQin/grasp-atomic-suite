!***********************************************************************
!                                                                      *
      subroutine lagcon(j, nprocs)
!                                                                      *
!   This  routine  includes  the Lagrange multiplier contribution in   *
!   the 'exchange' term.                                               *
!
!   Parameter nprocs is added so that it works for mpi and serial
!   programs.
!                                                                      *
!                                           Last update: 08 Dec 1992   *
!   Modified by Xinghong He                 Last update: 17 Aug 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: orbital_index_radix
      use def_c
      use grid_c
      use lagr_c
      use pote_c
      use scf_c
      use wave_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: j
      integer , intent(in) :: nprocs
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, iecck, l1, l2, m, mfm, i
      real(double) :: eps, wb, wa, wari
!-----------------------------------------------
!
!
      if (nec == 0) return
!
!ww      EPS = ACCY*0.1D 00
      eps = accy*0.01d00
!
!   Add contributions from off-diagonal parameters to exchange
!
      wb = 1.0d00/(ucf(j)*speed_of_light_au)/nprocs
      do k = 1, nec
!
!   Decode index
!
         iecck = iecc(k)
         l1 = iecck/key
         l2 = iecck - key*l1
!
         if (j == l1) then
            m = l2
         else if (j == l2) then
            m = l1
         else
            cycle
         endif
!
         wa = ecv(k)*wb
         if (abs(wa) < eps) cycle
!
!   ADD CONTRIBUTIONS TO EXCHANGE TERMS
!
         mfm = orbital_last_grid_index(m)
         do i = 1, mfm
            wari = wa*radial_grid(i)
            xp(i) = xp(i) + wari*radial_small_components(i,m)
            xq(i) = xq(i) - wari*radial_large_components(i,m)
         end do
!
      end do
!
      return
      end subroutine lagcon
