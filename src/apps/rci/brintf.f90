!***********************************************************************
!                                                                      *
      real(real64) function brintf (itype, orbital_a, orbital_b, orbital_c, orbital_d, k)
!                                                                      *
!   Computes integrals for the transverse photon interaction.          *
!                                                                      *
!   Call(s) to: [RCI92]: BESSEL, BRRA.                                 *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use grasp_parameters,   only: orbital_index_radix
      use breit_bessel_cache_keys, only: keep
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use bessel_i, only: bessel
      use brra_i, only: brra
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: itype
      integer, intent(in) :: orbital_a
      integer, intent(in) :: orbital_b
      integer, intent(in) :: orbital_c
      integer, intent(in) :: orbital_d
      integer, intent(in) :: k
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: orbital_label_radix = orbital_index_radix
      integer, parameter :: key2 = orbital_label_radix*orbital_label_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icod, icod1, icod2, icod3, icod4
!-----------------------------------------------
!
      goto (1,4,3,6,2,5) itype
!
!   Type 1 and 5 integrals require j(k), n(k) Bessel fuctions
!   Type 5 integrals only require w = wab Bessel functions
!
    1 continue
      if (orbital_a/=orbital_c .or. orbital_b/=orbital_d .or. orbital_a/=orbital_d .or. orbital_c==orbital_b) call bessel ( &
         & orbital_c, orbital_d, 1, 2&
         , k)
    2 continue
      call bessel (orbital_a, orbital_b, 1, 1, k)
      goto 6
!
!   Type 3 integrals require j(k), n(k) Bessel functions for either
!   w = wab or w = cd whichever is non-zero.
!
    3 continue
      if (orbital_a /= orbital_b) call bessel (orbital_a, orbital_b, 1, 1, k)
      if (orbital_c /= orbital_d) call bessel (orbital_c, orbital_d, 1, 2, k)
      goto 6
!
!   Type 2 and 6 integrals require j(k), n(k) and j(k+2), n(k+2)
!   Bessel fuctions
!   Type 6 integrals only require w = wab Bessel functions.
!
    4 continue
      if (orbital_a/=orbital_c .or. orbital_b/=orbital_d .or. orbital_a/=orbital_d .or. orbital_c/=orbital_b) then
!
         icod = max(orbital_c,orbital_d) + orbital_label_radix*min(orbital_c,orbital_d)
         icod1 = icod + key2*(k - 1)
         icod2 = icod + key2*(k + 1)
         icod = max(orbital_a,orbital_b) + orbital_label_radix*min(orbital_a,orbital_b)
         icod3 = icod + key2*(k - 1)
         icod4 = icod + key2*(k + 1)
         if (icod1==keep(1,2) .and. icod2==keep(2,2) .and. icod3==keep(1,1)&
             .and. icod4==keep(2,1)) goto 6
         if (icod1==keep(1,1) .and. icod2==keep(2,1) .and. icod3==keep(1,2)&
             .and. icod4==keep(2,2)) goto 6
         call bessel (orbital_c, orbital_d, 1, 2, k - 1)
         call bessel (orbital_c, orbital_d, 2, 2, k + 1)
      end if
!
    5 continue
      call bessel (orbital_a, orbital_b, 1, 1, k - 1)
      call bessel (orbital_a, orbital_b, 2, 1, k + 1)
!
!   Compute the integral
!
    6 continue
      brintf = brra(itype,orbital_a,orbital_b,orbital_c,orbital_d,k)
!
      return
      end function brintf
