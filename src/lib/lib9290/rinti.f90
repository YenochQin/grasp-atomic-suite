!***********************************************************************
!                                                                      *
      real(real64) function rinti (j, k, mode)
!                                                                      *
!   The value of this  function is the one-electron integral I (J,K)   *
!   for  orbitals  J, K. The analytical expression for this quantity   *
!   is given as  eq (9) in  I P Grant, B J McKenzie, P H Norrington,   *
!   D F Mayers, and N C Pyper,  Computer  Phys Commun 21 (1980) 211 .  *
!                                                                      *
!   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use debug_c
      use def_c,           only: speed_of_light_au
      use grid_c
      use npot_c,          only: negative_radius_times_potential
      use orb_c
      use tatb_c,          only: ta, tb, mtp
      use wave_c,          only: orbital_last_grid_index,radial_large_components,radial_small_components
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dpbdt_i
      use quad_i
      implicit none (type, external)
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: mode
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(real64) :: piece1, piece2, piece3, piece4
!-----------------------------------------------
!
!   Stop if orbitals J and K have different kappa values
!
      if (orbital_kappa(j) /= orbital_kappa(k)) then
         write (*, 300) orbital_principal_n(j), orbital_symmetry_labels(j), orbital_principal_n(k), orbital_symmetry_labels(k)
         stop
      endif
!
      mtp = max(orbital_last_grid_index(j),orbital_last_grid_index(k))
!
!   Kinetic energy contribution
!
!   Piece involving derivatives
!
      call dpbdt (k)
      ta(1) = 0._real64
      do i = 2, mtp
         ta(i) = radial_small_components(i,j)*ta(i) - radial_large_components(i,j)*tb(i)
      end do
      call quad (piece1)
      piece1 = speed_of_light_au*piece1/transformed_grid_step
!
!   Pieces not involving derivatives
!
      ta(1) = 0._real64
      do i = 2, mtp
         ta(i) = radial_grid_derivative(i)*radial_small_components(i,j)*radial_small_components(i,k)
      end do
      call quad (piece2)
      piece2 = -2._real64*speed_of_light_au*speed_of_light_au*piece2
!
      ta(1) = 0._real64
      do i = 2, mtp
         ta(i) = radial_derivative_over_radius(i)*(radial_large_components(i,j)*radial_small_components(i,k) + &
            & radial_small_components(i,j)*radial_large_components(i,k))
      end do
      call quad (piece3)
      piece3 = piece3*speed_of_light_au*dble(orbital_kappa(k))
!
!   Contribution from nuclear potential only if MODE is 0
!
      if (mode == 0) then
!
         ta(1) = 0._real64
         do i = 2, mtp
            ta(i) = radial_derivative_over_radius(i)*negative_radius_times_potential(i)*(radial_large_components(i,j)* &
               & radial_large_components(i,k) + radial_small_components(i,j)*radial_small_components(i,k))
         end do
         call quad (piece4)
         piece4 = -piece4
!
      else
         piece4 = 0._real64
      endif
!
      rinti = piece1 + piece2 + piece3 + piece4
!
!   Debug printout
!
      if (mode==0 .and. ldbpr(4)) write (99, 301) orbital_principal_n(j), orbital_symmetry_labels(j), &
         & orbital_principal_n(k), orbital_symmetry_labels(k), &
         rinti
      if (mode/=0 .and. ldbpr(5)) write (99, 302) orbital_principal_n(j), orbital_symmetry_labels(j), &
         & orbital_principal_n(k), orbital_symmetry_labels(k), &
         rinti
!
      return
!
  300 format('RINTI: Attempt to calculate I(',1i2,1a2,',',1i2,1a2,')')
  301 format(/,' I (',1i2,1a2,',',1i2,1a2,') = ',1p,d19.12,/)
  302 format(/,' K (',1i2,1a2,',',1i2,1a2,') = ',1p,d19.12,/)
      return
!
      end function rinti
