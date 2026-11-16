!***********************************************************************
!                                                                      *
      real(real64) function vinti (j, k)
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
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
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
      integer, intent(in) :: k
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, kpj, kpk, ifact1, ifact2
      real(real64) :: piece1, fact1, fact2, piece2
!-----------------------------------------------
!
      mtp = max(orbital_last_grid_index(j),orbital_last_grid_index(k))
!
!   Piece involving derivatives
!
      call dpbdt (k)
      ta(1) = 0.0_real64
      do i = 2, mtp
         ta(i) = radial_large_components(i,j)*ta(i) + radial_small_components(i,j)*tb(i)
      end do
      call quad (piece1)
      piece1 = piece1/transformed_grid_step
!
!   Pieces not involving derivatives
!
      kpj = orbital_kappa(j)
      kpk = orbital_kappa(k)
      ifact1 = kpj*(kpj + 1) - kpk*(kpk + 1)
      fact1 = 0.5_real64*dble(ifact1)
      ifact2 = (-kpj*((-kpj) + 1)) + kpk*((-kpk) + 1)
      fact2 = 0.5_real64*dble(ifact2)
      ta(1) = 0.0_real64
      do i = 2, mtp
         ta(i) = radial_derivative_over_radius(i)*(fact1*radial_large_components(i,j)*radial_large_components(i,k) + &
            & fact2*radial_small_components(i,j)*radial_small_components(i,k))
      end do
      call quad (piece2)
!
      vinti = piece1 - piece2
!
!   Debug printout
!
      if (ldbpr(6)) write (99, 300) orbital_principal_n(j), orbital_symmetry_labels(j), orbital_principal_n(k), &
         & orbital_symmetry_labels(k), vinti
!
      return
!
  300 format(/,'VINTI: V (',1i2,1a2,',',1i2,1a2,') = ',1p,d19.12)
      return
!
      end function vinti
