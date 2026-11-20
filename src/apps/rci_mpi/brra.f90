!***********************************************************************
!                                                                      *
      real(real64) function brra (itype, orbital_a, orbital_c, orbital_b, orbital_d, k)
!                                                                      *
!   This routine evaluates the transverse interaction integrals:       *
!                                                                      *
!      ITYPE = 1: General R (k; a c | b d )                            *
!            = 2: General S (k; a c | b d )                            *
!            = 3: R (k; a , b d )                                      *
!            = 4: F (k; a , b )                                        *
!            = 5: G (k; a , b ) type integral                          *
!            = 6: H (k; a , b ) type integral                          *
!                                                                      *
!   Call(s) to: [RCI92]: RKINT, SKINT.                                 *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_radial_points
      use numeric_constants, only: eleven, eps, four, half, one, seven, ten, tenth, three, two, zero
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rkint_i, only: rkint
      use skint_i, only: skint
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: itype
      integer, intent(in) :: orbital_a
      integer, intent(in) :: orbital_c
      integer, intent(in) :: orbital_b
      integer, intent(in) :: orbital_d
      integer  :: k
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mxrac, i, mxrbd
      real(real64) :: rac(max_radial_points), rbd(max_radial_points)
!-----------------------------------------------
!
      mxrac = min(orbital_last_grid_index(orbital_a),orbital_last_grid_index(orbital_c))
      do i = 1, mxrac
         rac(i) = radial_large_components(i,orbital_a)*radial_small_components(i,orbital_c)
      end do
!
      mxrbd = min(orbital_last_grid_index(orbital_b),orbital_last_grid_index(orbital_d))
      do i = 1, mxrbd
         rbd(i) = radial_large_components(i,orbital_b)*radial_small_components(i,orbital_d)
      end do
!
      goto (21,22,23,24,25,26) itype
!
!   ITYPE = 1
!
   21 continue
      if (orbital_a==orbital_b .and. orbital_c==orbital_d) goto 9
      if (orbital_a==orbital_d .and. orbital_c==orbital_b) goto 10
      brra = (rkint(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,1) + rkint(rac,orbital_a,orbital_c,rbd,orbital_b, &
         & orbital_d,k,2)&
          + rkint(rbd,orbital_b,orbital_d,rac,orbital_a,orbital_c,k,1) + rkint(rbd,orbital_b,orbital_d,rac,orbital_a,orbital_c, &
             & k,2))*&
         half
      return
!
!   ITYPE = 2
!
   22 continue
      if (orbital_a==orbital_b .and. orbital_c==orbital_d) goto 26
      if (orbital_a==orbital_d .and. orbital_c==orbital_b) goto 26
      brra = (skint(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,1) + skint(rac,orbital_a,orbital_c,rbd,orbital_b, &
         & orbital_d,k,2))*&
         half
      return
!
!   ITYPE = 3
!
   23 continue
      if (orbital_a == orbital_c) then
         do i = 1, mxrbd
            rbd(i) = rbd(i) + radial_large_components(i,orbital_d)*radial_small_components(i,orbital_b)
         end do
         brra = (rkint(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,0) + rkint(rbd,orbital_b,orbital_d,rac,orbital_a, &
            & orbital_c,k,0&
            ) + rkint(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,2) + rkint(rbd,orbital_b,orbital_d,rac,orbital_a, &
               & orbital_c,k,2)&
            )*half
         return
      end if
      do i = 1, mxrac
         rac(i) = rac(i) + radial_large_components(i,orbital_c)*radial_small_components(i,orbital_a)
      end do
      brra = (rkint(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,1) + rkint(rbd,orbital_b,orbital_d,rac,orbital_a, &
         & orbital_c,k,1)&
          + rkint(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,0) + rkint(rbd,orbital_b,orbital_d,rac,orbital_a,orbital_c, &
             & k,0))*&
         half
      return
!
!   ITYPE = 4
!
   24 continue
      brra = rkint(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,0) + rkint(rbd,orbital_b,orbital_d,rac,orbital_a,orbital_c,k,0)
      return
!
!   ITYPE = 5
!
   25 continue
      if (orbital_a==orbital_d .and. orbital_c==orbital_b) goto 10
    9 continue
      brra = two*rkint(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,1)
      return
   10 continue
      brra = rkint(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,1) + rkint(rbd,orbital_b,orbital_d,rac,orbital_a,orbital_c,k,1)
      return
!
!   ITYPE = 6
!
   26 continue
      brra = skint(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,1)
      return
!
      end function brra
